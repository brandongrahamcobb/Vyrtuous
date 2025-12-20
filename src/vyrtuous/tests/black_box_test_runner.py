''' black_box_test_runner.py The purpose of this program is to provide the black box test runner.
    Copyright (C) 2025  https://gitlab.com/vyrtuous/vyrtuous

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''
from vyrtuous.tests.test_admin_commands import AdminCommandsTest
# from vyrtuous.tests.aliases_test import AliasesTest
# from vyrtuous.tests.coord_commands_test import CoordinatorCommandsTest
# from vyrtuous.tests.dev_commands_test import DeveloperCommandsTest
# from vyrtuous.tests.events_test import EventsTest
# from vyrtuous.tests.everyone_commands_test import EveryoneCommandsTest
# from vyrtuous.tests.help_command_test import HelpCommandTest
# from vyrtuous.tests.mod_commands_test import ModeratorCommandsTest
# from vyrtuous.tests.owner_commands_test import OwnerCommandsTest
# from vyrtuous.tests.tasks_test import TasksTest

class BlackBoxTestRunner():

    @classmethod
    async def run(cls, client):
        await AdminCommandsTest(client).begin()
        # await AliasesTest.begin()
        # await CoordinatorCommandsTest.begin()
        # await DeveloperCommandsTest.begin()
        # await EventsTest.begin()
        # await EveryoneCommandsTest.begin()
        # await HelpCommandTest.begin()
        # await ModeratorCommandsTest.begin()
        # await OwnerCommandsTest.begin()
        # await TasksTest.begin()