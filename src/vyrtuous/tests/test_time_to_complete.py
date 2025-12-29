
# ''' test_time_to_complete.py The purpose of this program is to provide the tests for the TimeToComplete module.
#     Copyright (C) 2025  https://gitlab.com/vyrtuous/vyrtuous

#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.

#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.

#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.
# '''
# from typing import Optional
# from vyrtuous.utils.time_to_complete import TimeToComplete
# import pytest
# import time

# @pytest.fixture
# def end():
#     return time.perf_counter()

# @pytest.fixture
# def start():
#     return time.perf_counter()

# def test_time_elapsed_measurement(start: Optional[float], end: Optional[float]):
#     try:
#         is_around_one_second = 1.0
#         counter = TimeToComplete()
#         elapsed = counter.time_elapsed_measurement(start, end)
#         if time.perf_counter() == start:
#             start = time.perf_counter()
#         if start > end:
#             is_around_one_second = abs(start - end)
#         else:
#             is_around_one_second = start - end
#         if not counter.is_around_one_second(is_around_one_second):
#             is_around_one_second = 1
#             time.sleep(1)
#             end = time.perf_counter()
#             elapsed = counter.time_elapsed_measurement(start, end)
#             is_around_one_second = elapsed
#         assert counter.is_around_one_second(is_around_one_second)
#     except TypeError:
#         raise TypeError("Fixtures end and start must be floats representing time in seconds.")
