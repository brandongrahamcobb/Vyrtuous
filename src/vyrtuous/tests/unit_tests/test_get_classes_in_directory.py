from pathlib import Path

import pytest

from vyrtuous.inc.helpers import DIR_BASE
from vyrtuous.utils.dir_to_classes import dir_to_classes


def test_get_clss_in_dir():
    dbd = Path(__file__).resolve().parents[2] / "database"
    clss = dir_to_classes(str(dbd))
    clss_str = "\n".join([str(cls) for cls in clss])
    print(clss_str)
