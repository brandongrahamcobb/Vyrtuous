from vyrtuous.fields.duration import DurationObject


def update_duration_display(td):
    if td is None:
        label = "Never"
    total_seconds = int(td.total_seconds())
    days, remainder = divmod(total_seconds, 86400)
    hours, remainder = divmod(remainder, 3600)
    minutes, _ = divmod(remainder, 60)
    parts = []
    if days:
        parts.append(f"in {days}d")
    if hours:
        parts.append(f"in {hours}h")
    if minutes:
        parts.append(f"in {minutes}m")
    label = " ".join(parts)
    return label


if __name__ == "__main__":
    td = DurationObject("+8h").to_timedelta()
    print(update_duration_display(td))


class AliasModal(discord.ui.Modal):

    def __init__(self, partial_action_information, state):
        super().__init__(
            title=f'{partial_action_information["alias_class"].SINGULAR} Reason'
        )
        self.action_information = partial_action_information
        self.bot = DiscordBot.get_instance()
        existing_reason = (
            self.action_information.get("action_existing").reason
            if self.action_information.get("action_existing")
            else None
        )
        self.reason = discord.ui.TextInput(
            label=f'Type {self.action_information["alias_class"].SINGULAR.lower()} reason',
            style=discord.TextStyle.paragraph,
            required=True,
            default=existing_reason or "",
        )
        self.add_item(self.reason)
        self.state = state

    async def on_submit(self, interaction):
        try:
            await has_equal_or_lower_role_wrapper(
                source=interaction,
                member_snowflake=self.action_information["action_member_snowflake"],
                sender_snowflake=interaction.user.id,
            )
        except HasEqualOrLowerRole as e:
            return await interaction.response.send_message(str(e), ephemeral=True)
        executor_role = await has_equal_or_lower_role_wrapper(
            source=interaction,
            member_snowflake=self.action_information["action_member_snowflake"],
            sender_snowflake=interaction.user.id,
        )
        self.action_information["action_executor_role"] = (
            executor_role  # Adds executor role
        )
        expires_in_timedelta = self.action_information["action_duration"].to_timedelta()
        if self.action_information["action_duration"].number != 0:
            if expires_in_timedelta.total_seconds() < 0:
                return await self.state.end(
                    warning="You are not authorized to decrease "
                    "the duration below the current time."
                )
        if (
            "action_existing" in self.action_information.keys()
            and self.action_information["action_existing"]
        ):
            if (
                expires_in_timedelta.total_seconds()
                > self.action_information["action_channel_cap"]
                or self.action_information["action_duration"].number == 0
            ):
                if executor_role == "Moderator":
                    duration_str = DurationObject.from_seconds(
                        self.action_information["action_channel_cap"]
                    )
                    return await self.state.end(
                        warning=f"Cannot set the {self.action_information["alias_class"].SINGULAR} beyond {duration_str} as a "
                        f"{executor_role} in {channel_obj.mention}."
                    )
        self.action_information["action_reason"] = self.reason.value
        channel_obj = interaction.guild.get_channel(
            self.action_information.get("action_channel_snowflake", None)
        )
        member_obj = interaction.guild.get_member(
            self.action_information["action_member_snowflake"]
        )
        func = self.action_information["alias_class"].get_handler()
        from types import SimpleNamespace

        alias = SimpleNamespace(
            alias_name=None,
            alias_type="ban",
            channel_snowflake=self.action_information.get("action_channel_snowflake", None),
        )
        await func(
            alias=alias,
            action_information=self.action_information,
            channel=channel_obj,
            member=member_obj,
            message=interaction.message,
            state=self.state,
        )
        return await interaction.response.defer()


class ChannelView(discord.ui.View):
    def __init__(self, partial_action_information, interaction, state):
        super().__init__(timeout=120)
        self.action_information = {}
        self.action_information.update(
            partial_action_information
        )  # Adds alias_class, action_guild_snowflake, action_member_snowflake
        self.duration_buttons_added = False
        self.interaction = interaction
        self.channel_select.options = self._build_channel_options()
        self.expires_in_timedelta = None
        self.state = state

    def _build_channel_options(self):
        guild = self.interaction.guild
        if not guild:
            return []
        return [
            discord.SelectOption(label=ch.name, value=str(ch.id))
            for ch in guild.channels
            if isinstance(ch, discord.VoiceChannel)
        ]

    @discord.ui.select(
        placeholder="Select channel",
        options=[],
    )
    async def channel_select(self, interaction, select):
        channel = interaction.guild.get_channel(int(select.values[0]))
        self.action_information.get("action_channel_snowflake", None) = (
            channel.id
        )  # Adds channel_snowflake
        action_channel_cap = await Alias.generate_cap_duration(
            channel_snowflake=channel.id,
            guild_snowflake=interaction.guild.id,
            moderation_type=self.action_information["alias_class"].ACT,
        )
        self.action_information["action_channel_cap"] = (
            action_channel_cap  # Adds channel cap
        )
        self.channel_select.placeholder = channel.name
        await interaction.response.defer()
        where_kwargs = {
            "channel_snowflake": self.action_information.get("action_channel_snowflake", None),
            "guild_snowflake": interaction.guild.id,
            "member_snowflake": self.action_information["action_member_snowflake"],
        }
        action_existing = await self.action_information["alias_class"].select(
            **where_kwargs, singular=True
        )
        self.action_information["action_modification"] = False
        if action_existing:
            self.action_information["action_existing"] = action_existing
            self.action_information["action_modification"] = True
            self.update_duration_display(action_existing.expires_in)
        if not self.duration_buttons_added:
            self._add_duration_buttons()
            self.duration_buttons_added = True
        await interaction.edit_original_response(view=self)

    def _add_duration_buttons(self):
        self.add_item(DurationButton("Reason Only"))
        self.add_item(DurationButton("+1h"))
        self.add_item(DurationButton("+8h"))
        self.add_item(DurationButton("+24h"))
        self.add_item(DurationButton("Permanent"))

    async def on_timeout(self):
        for item in self.children:
            item.disabled = True

    def update_duration_display(self, dt):
        td = DurationObject.from_expires_in(dt).to_timedelta()
        if td is None:
            label = "Never"
        total_seconds = int(td.total_seconds())
        days, remainder = divmod(total_seconds, 86400)
        hours, remainder = divmod(remainder, 3600)
        minutes, _ = divmod(remainder, 60)
        parts = []
        if days:
            parts.append(f"in {days}d")
        if hours:
            parts.append(f"in {hours}h")
        if minutes:
            parts.append(f"in {minutes}m")
        label = " ".join(parts)
        self.add_item(DurationDisplayButton(label))


class DurationButton(discord.ui.Button):
    def __init__(self, value):
        if value.startswith("+"):
            style = discord.ButtonStyle.primary
        elif value == "Permanent":
            style = discord.ButtonStyle.danger
        else:
            style = discord.ButtonStyle.secondary
        super().__init__(label=value, style=style)
        self.value = value

    async def callback(self, interaction):
        if not self.view.action_information.get("action_channel_snowflake"):
            return await interaction.response.send_message(
                "Select an option first.", ephemeral=True
            )
        if self.value != "Reason Only":
            if self.value == "Permanent":
                self.value = 0
            action_duration = DurationObject(self.value)
            self.view.action_information["action_duration"] = (
                action_duration  # Adds duration
            )
            self.view.expires_in_timedelta = action_duration.to_timedelta()
            if (
                "action_existing" in self.view.action_information.keys()
                and self.view.action_information["action_existing"]
            ):
                action_expires_in = (
                    self.view.action_information["action_existing"].expires_in
                    + self.view.expires_in_timedelta
                )
            else:
                action_expires_in = (
                    datetime.now(timezone.utc) + self.view.expires_in_timedelta
                )
            self.view.action_information["action_expires_in"] = action_expires_in
            action_expires_in = (
                datetime.now(timezone.utc) + self.view.expires_in_timedelta
            )
            self.view.action_information["action_expires_in"] = (
                action_expires_in  # Adds expires in
            )
        modal = AliasModal(
            partial_action_information=self.view.action_information,
            state=self.view.state,
        )
        await interaction.response.send_modal(modal)


class DurationDisplayButton(discord.ui.Button):
    def __init__(self, label):
        super().__init__(
            label=label, style=discord.ButtonStyle.secondary, disabled=True
        )

    @app_commands.command(name="ban", description="Set a ban.")
    @app_commands.describe(
        member="Specify member ID or mention.",
    )
    @moderator_predicator()
    async def ban_app_command(
        self, interaction: discord.Interaction, member: AppMemberSnowflake
    ):
        state = StateService(source=interaction)
        do = DiscordObject(interaction=interaction)
        member_dict = await do.determine_from_target(target=member)
        partial_action_information = {
            "alias_class": Ban,
            "action_guild_snowflake": interaction.guild.id,
            "action_member_snowflake": member_dict["id"],
        }
        view = ChannelView(
            partial_action_information=partial_action_information,
            interaction=interaction,
            state=state,
        )
        return await interaction.response.send_message(
            content=f"Banning {member_dict['mention']}", view=view
        )
