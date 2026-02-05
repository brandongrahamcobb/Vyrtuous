"""dir_to_classes.py The purpose of this program is to provide the dir_to_classes utility module.

Copyright (C) 2025  https://github.com/brandongrahamcobb/Vyrtuous.git

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

import importlib.util
import inspect
import sys
from pathlib import Path

from vyrtuous.utils.logger import logger


def dir_to_classes(dir_paths, parent):
    classes = []
    for dir_path in dir_paths:
        for py_file in dir_path.rglob("*.py"):
            module_name = py_file.stem
            spec = importlib.util.spec_from_file_location(module_name, str(py_file))
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            for _, cls in inspect.getmembers(module, inspect.isclass):
                if cls.__module__ != module.__name__:
                    continue
                if getattr(cls, "__skip_db_discovery__", False):
                    continue
                if issubclass(cls, parent) and cls is not parent:
                    classes.append(cls)
    return classes


def skip_db_discovery(cls):
    cls.__skip_db_discovery__ = True
    return cls
