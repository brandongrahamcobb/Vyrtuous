import ast
import os


def collect_method_names():
    method_names = set()
    paths = [p for p in os.listdir(".") if p.endswith(".py")]
    for path in paths:
        with open(path, "r", encoding="utf-8") as file:
            try:
                tree = ast.parse(file.read())
            except SyntaxError:
                continue
        for node in ast.walk(tree):
            if isinstance(node, ast.FunctionDef):
                method_names.add(node.name)
            if isinstance(node, ast.AsyncFunctionDef):
                method_names.add(node.name)
    return method_names


def main():
    method_names = collect_method_names()
    for name in sorted(method_names):
        print(name)


main()
