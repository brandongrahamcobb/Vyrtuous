from copy import copy

import discord
from discord.ext import commands

from vyrtuous.upload.upload import Upload


class UploadService:
    MODEL = Upload

    def __init__(
        self,
        *,
        bot,
        database_factory=None,
    ):
        self.__bot = bot
        self.__source = None
        self.__database_factory = copy(database_factory)
        self.__database_factory.model = self.MODEL

    def __extract_command_and_args(self):
        if isinstance(self.__source, discord.Interaction):
            command_name = getattr(
                getattr(self.__source, "command", None), "name", None
            )
            namespace = getattr(self.__source, "namespace", None)
            arguments = ""
            if namespace:
                arguments = " ".join(
                    str(value)
                    for value in vars(namespace).values()
                    if value is not None
                )
            return command_name, arguments
        elif isinstance(self.__source, commands.Context):
            command_name = self.__source.command.name if self.__source.command else None
            arguments = " ".join(map(str, self.__source.args[1:]))
            return command_name, arguments
        elif isinstance(self.__source, discord.Message):
            prefix = self.__bot.config["discord_command_prefix"]
            content = self.__source.content.strip()
            if content.startswith(prefix):
                content = content[len(prefix) :]
            parts = content.split()
            if not parts:
                return None, None
            command_name = parts[0]
            arguments = " ".join(parts[1:]) if len(parts) > 1 else ""
            return command_name, arguments
        return None, None

    async def __send_prompt(self):
        if isinstance(self.__source, discord.Interaction):
            await self.__source.followup.send(
                "Upload a file in your next message.",
                ephemeral=False,
            )
            return self.__source.channel, self.__source.user
        elif isinstance(self.__source, commands.Context):
            await self.__source.send("Upload a file in your next message.")
            return self.__source.channel, self.__source.author
        elif isinstance(self.__source, discord.Message):
            await self.__source.channel.send("Upload a file in your next message.")
            return self.__source.channel, self.__source.author
        return None, None

    async def __wait_for_upload(self, channel, user):
        def check(m: discord.Message):
            return m.channel == channel and m.author == user and len(m.attachments) > 0

        return await self.__bot.wait_for("message", timeout=300, check=check)

    async def __save_upload(self, attachment, command_name, tag):
        file_bytes = await attachment.read()
        obj = Upload(
            command_name=command_name,
            file_bytes=file_bytes,
            filename=attachment.filename,
            tag=tag,
        )
        await self.__database_factory.upsert(obj)

    async def request_upload(self, source):
        self.__source = source
        if not self.__source:
            return False
        command_name, arguments = self.__extract_command_and_args()
        if not command_name:
            return False
        channel, user = await self.__send_prompt()
        if not channel:
            return False
        try:
            message = await self.__wait_for_upload(channel, user)
        except Exception:
            return False
        attachment = message.attachments[0]
        tag = message.content
        try:
            await self.__save_upload(attachment, command_name, tag)
        except Exception:
            return False
        return True

    async def build_latex_document(self):
        import re

        def escape_latex(text):
            return re.sub(r"([&_#%$])", r"\\\1", str(text))

        import os

        os.makedirs("images", exist_ok=True)
        async with self.__bot.db_pool.acquire() as conn:
            rows = await conn.fetch(
                "select command_name,file_bytes,filename,tag,created_at,updated_at from uploads order by command_name,created_at"
            )
        with open("uploads.tex", "w") as f:
            f.write(
                r"\documentclass{article}\usepackage{graphicx}\usepackage{float}\usepackage{geometry}\usepackage{hyperref}\geometry{margin=1in}\title{Uploads Archive}\date{\today}\begin{document}\maketitle\tableofcontents\newpage"
            )
            current_command = None
            for row in rows:
                command_name = row["command_name"]
                tag = row["tag"]
                filename = row["filename"]
                created_at = row["created_at"]
                updated_at = row["updated_at"]
                file_bytes = row["file_bytes"]
                if command_name != current_command:
                    if current_command is not None:
                        f.write(r"\clearpage")
                    f.write(rf"\section{{Command: {escape_latex(command_name)}}}")
                    current_command = command_name
                f.write(rf"\subsection{{Tag: {escape_latex(tag)}}}")
                f.write(rf"\textbf{{Filename:}} {escape_latex(filename)}\\")
                f.write(rf"\textbf{{Created:}} {created_at}\\")
                f.write(rf"\textbf{{Updated:}} {updated_at}")
                image_path = f"images/{command_name}_{tag}_{filename}"
                image_path = re.sub(r"\s+", "_", image_path)
                with open(image_path, "wb") as img:
                    img.write(file_bytes)
                f.write(
                    rf"\begin{{figure}}[H]\centering\includegraphics[width=\linewidth,height=0.50\textheight,keepaspectratio]{{{image_path}}}\end{{figure}}"
                )
            f.write(r"\end{document}")
