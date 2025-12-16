import pytest
import time

def test_long_task():
    assert time.sleep(5)

def test_sleep_broken():
    assert time.sleep(-1)

def test_short_task(*args):
    assert time.sleep(5)