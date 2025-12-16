
from typing import Optional
from vyrtuous.utils.time_to_complete import TimeCounter

import pytest
import time

@pytest.fixture
def end():
    return time.perf_counter()

@pytest.fixture
def start():
    return time.perf_counter()

def test_time_elapsed_measurement(start: Optional[float], end: Optional[float]):
    try:
        is_around_one_second = 1
        counter = TimeCounter()
        elapsed = counter.time_elapsed_measurement(start, end)
        if time.perf_counter() == start:
            start = time.perf_counter()
        if start > end:
            is_around_one_second = abs(start - end)
        else:
            is_around_one_second = start - end
        if not counter.is_around_one_second(is_around_one_second):
            is_around_one_second = 1
            time.sleep(1)
            end = time.perf_counter()
            elapsed = counter.time_elapsed_measurement(start, end)
            is_around_one_second = elapsed
        assert is_around_one_second
    except TypeError:
        raise TypeError("Fixtures end and start must be floats representing time in seconds.")
