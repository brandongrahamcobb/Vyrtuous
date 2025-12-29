import os
import fileinput
import sys

if len(sys.argv) < 2:
    print("Usage: python replace_imports.py <folder_path> [--recursive]")
    sys.exit(1)

folder_path = sys.argv[1]
recursive = "--recursive" in sys.argv

def process_file(filepath):
    for line in fileinput.input(filepath, inplace=True):
        if line.startswith("from vyrtuous.tests.black_box.test_suite import"):
            print("from vyrtuous.tests.black_box.test_suite import *")
        else:
            print(line, end="")

if recursive:
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith(".py"):
                process_file(os.path.join(root, file))
else:
    for file in os.listdir(folder_path):
        if file.endswith(".py"):
            process_file(os.path.join(folder_path, file))
