import pytest

def test_count_args(*args):
    assert len(args) <= 3