from pathlib import Path
import importlib.util
import inspect

def db_to_classes():
    classes = []
    dir_paths = []
    dir_paths.append(Path(__file__).resolve().parents[1] / "db/actions")
    dir_paths.append(Path(__file__).resolve().parents[1] / "db/mgmt")
    dir_paths.append(Path(__file__).resolve().parents[1] / "db/roles")
    dir_paths.append(Path(__file__).resolve().parents[1] / "db/rooms")
    for dir_path in dir_paths:
        for py_file in dir_path.rglob("*.py"):
            module_name = py_file.stem
            spec = importlib.util.spec_from_file_location(module_name, str(py_file))
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            for _, cls in inspect.getmembers(module, inspect.isclass):
                if cls.__module__ != module.__name__:
                    continue
                if getattr(cls, '__skip_db_discovery__', False):
                    continue
                classes.append(cls)
    return classes

def skip_db_discovery(cls):
    cls.__skip_db_discovery__ = True
    return cls
