from pathlib import Path

import pytest


from vyrtuous.inc.helpers import DIR_BASE
from vyrtuous.utils.get_clss_in_dir import get_clss_in_dir


def test_get_clss_in_dir():
    dbd = Path(__file__).resolve().parents[2] / "database"
    clss = get_clss_in_dir(str(dbd))
    clss_str = "\n".join([str(cls) for cls in clss])
    print(clss_str)
