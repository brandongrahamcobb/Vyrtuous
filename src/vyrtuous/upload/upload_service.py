import discord
from discord.ext import commands

from vyrtuous.upload.upload import Upload


class UploadService:
    def __init__(
        self,
        *,
        bot,
        database_factory=None,
    ):
        self.__bot = bot
        self.__source = None
        self.__database_factory = database_factory

    def __extract_command_and_args(self):
        if isinstance(self.__source, discord.Interaction):
            command_name = self.__source.command.name if self.__source.command else None
            arguments = " ".join(
                str(value)
                for value in getattr(self.__source, "namespace", {}).values()
                if value is not None
            )
            return command_name, arguments
        if isinstance(self.__source, commands.Context):
            command_name = self.__source.command.name if self.__source.command else None
            arguments = " ".join(
                map(str, self.__source.args[1:])
            )  # skip command itself
            return command_name, arguments
        if isinstance(self.__source, discord.Message):
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
            await self.__source.response.send_message(
                "Upload a file in your next message.",
                ephemeral=True,
            )
            return self.__source.channel, self.__source.user
        elif isinstance(self.__source, discord.Message) or isinstance(
            self.__source, commands.Context
        ):
            await self.__source.send("Upload a file in your next message.")
            return self.__source.channel, self.__source.author
        return None, None

    async def __wait_for_upload(self, channel, user):
        def check(m: discord.Message):
            return m.channel == channel and m.author == user and len(m.attachments) > 0

        return await self.__bot.wait_for("message", timeout=300, check=check)

    async def __save_upload(self, attachment, command_name, arguments):
        file_bytes = await attachment.read()
        obj = Upload(
            command_name=command_name,
            arguments=arguments,
            file_bytes=file_bytes,
            filename=attachment.filename,
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
        try:
            await self.__save_upload(attachment, command_name, arguments)
        except Exception:
            return False
        return True
