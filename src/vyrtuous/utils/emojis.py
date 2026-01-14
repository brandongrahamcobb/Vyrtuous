"""emojis.py A utility class for managing emojis used by the Vyrtuous bot.
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
"""

import random
EMOJIS = [
    "\U0001f436",
    "\U0001f431",
    "\U0001f42d",
    "\U0001f439",
    "\U0001f430",
    "\U0001f98a",
    "\U0001f43b",
    "\U0001f43c",
    "\U0001f428",
    "\U0001f42f",
    "\U0001f981",
    "\U0001f42e",
    "\U0001f437",
    "\U0001f43d",
    "\U0001f438",
    "\U0001f435",
    "\U0001f412",
    "\U0001f98d",
    "\U0001f9a7",
    "\U0001f414",
    "\U0001f427",
    "\U0001f426",
    "\U0001f424",
    "\U0001f423",
    "\U0001f425",
    "\U0001f986",
    "\U0001f9a2",
    "\U0001f989",
    "\U0001f99a",
    "\U0001f99c",
    "\U0001f43a",
    "\U0001f99d",
    "\U0001f9a8",
    "\U0001f9a1",
    "\U0001f417",
    "\U0001f434",
    "\U0001f984",
    "\U0001f41d",
    "\U0001f41b",
    "\U0001f98b",
    "\U0001f40c",
    "\U0001f41e",
    "\U0001f40c",
    "\U0001fab2",
    "\U0001f997",
    "\U0001f577",
    "\U0001f982",
    "\U0001f422",
    "\U0001f40d",
    "\U0001f98e",
    "\U0001f996",
    "\U0001f995",
    "\U0001f419",
    "\U0001f991",
    "\U0001f990",
    "\U0001f99e",
    "\U0001f980",
    "\U0001f421",
    "\U0001f420",
    "\U0001f41f",
    "\U0001f42c",
    "\U0001f988",
    "\U0001f433",
    "\U0001f40b",
    "\U0001f9ad",
    "\U0001f40a",
    "\U0001f406",
    "\U0001f405",
    "\U0001f403",
    "\U0001f402",
    "\U0001f42b",
    "\U0001f42a",
    "\U0001f999",
    "\U0001f992",
    "\U0001f98f",
    "\U0001f99b",
    "\U0001f418",
    "\U0001f998",
    "\U0001f9a5",
    "\U0001f9a6",
    "\U0001f9a8",
    "\U0001f9a9",
    "\U0001f54a",
]

def get_random_emoji(self):
    return random.choice(self.EMOJIS)
