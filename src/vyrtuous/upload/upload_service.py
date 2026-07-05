import discord
from discord.ext import commands

from vyrtuous.bot.discord_bot import DiscordBot
from vyrtuous.db.database_factory import DatabaseFactory
from vyrtuous.upload.upload import Upload

MODEL = Upload


def extract_command_and_args(source) -> tuple[str | None, str | None]:
    bot: DiscordBot = DiscordBot.get_instance()
    if isinstance(source, discord.Interaction):
        command_name = getattr(getattr(source, "command", None), "name", None)
        namespace = getattr(source, "namespace", None)
        arguments = ""
        if namespace:
            arguments = " ".join(
                str(value) for value in vars(namespace).values() if value is not None
            )
        return command_name, arguments
    elif isinstance(source, commands.Context):
        command_name = source.command.name if source.command else None
        arguments = " ".join(map(str, source.args[1:]))
        return command_name, arguments
    elif isinstance(source, discord.Message):
        prefix = bot.config["discord_command_prefix"]
        content = source.content.strip()
        if content.startswith(prefix):
            content = content[len(prefix) :]
        parts = content.split()
        if not parts:
            return None, None
        command_name = parts[0]
        arguments = " ".join(parts[1:]) if len(parts) > 1 else ""
        return command_name, arguments
    return None, None


async def send_prompt(
    source,
) -> tuple[
    discord.Message | None,
    discord.Member | discord.User | None,
]:
    if isinstance(source, discord.Interaction):
        await source.followup.send(
            "Upload a file in your next message.",
            ephemeral=False,
        )
        message: discord.Message | discord.interactions.InteractionMessage = (
            await source.original_response()
        )
        return message, source.user
    elif isinstance(source, commands.Context):
        message = await source.send("Upload a file in your next message.")
        return message, source.author
    elif isinstance(source, discord.Message):
        message = await source.channel.send("Upload a file in your next message.")
        return message, source.author
    return None, None


async def wait_for_upload(channel, user) -> discord.Message:
    bot: DiscordBot = DiscordBot.get_instance()

    def check(m: discord.Message):
        return m.channel == channel and m.author == user and len(m.attachments) > 0

    return await bot.wait_for("message", timeout=300, check=check)


async def save_upload(attachment, command_name, tag) -> None:
    database_factory: DatabaseFactory = DatabaseFactory(MODEL)
    file_bytes = await attachment.read()
    obj = Upload(
        command_name=command_name,
        file_bytes=file_bytes,
        filename=attachment.filename,
        tag=tag,
    )
    await database_factory.upsert(obj)


async def request_upload(source) -> bool:
    if not source:
        return False
    command_name, arguments = extract_command_and_args(source)
    if not command_name:
        return False
    message, user = await send_prompt(source)
    if message is None:
        return False
    try:
        message = await wait_for_upload(message.channel, user)
    except Exception:
        return False
    attachment = message.attachments[0]
    tag = message.content
    try:
        await save_upload(attachment, command_name, tag)
    except Exception:
        return False
    return True


async def build_latex_document() -> None:
    bot: DiscordBot = DiscordBot.get_instance()
    import re

    def escape_latex(text):
        return re.sub(r"([&_#%$])", r"\\\1", str(text))

    import os

    os.makedirs("images", exist_ok=True)
    async with bot.db_pool.acquire() as conn:
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
