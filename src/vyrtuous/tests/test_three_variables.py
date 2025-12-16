from vyrtuous.inc.helpers import *
import inspect
import pytest
import vyrtuous

def test_count_args(*args):
    assert len(args) <= 3

def test_verify_less_than_three_args():
    for name, obj in inspect.getmembers(vyrtuous):
        if inspect.isfunction(obj) or inspect.ismethod(obj):
            params = inspect.signature(obj).parameters
            assert len(params) <= 3, f"Function '{name}' has more than 3 parameters: {list(params.keys())}"