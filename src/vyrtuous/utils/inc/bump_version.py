import os
import re

def get_current_version(path):
    try:
        with open(path, "r") as f:
            content = f.read()
        match = re.search(r'export VYRTUOUS_VERSION="(\d+)\.(\d+)\.(\d+)"', content)
        if match:
            return tuple(map(int, match.groups())), content
    except FileNotFoundError:
        pass
    return (0, 0, 0), ""

def increment_version(version):
    major, minor, patch = version
    patch += 1
    if patch > 9:
        minor += 1
        patch = 0
    if minor > 9:
        major += 1
        minor = 0
    return major, minor, patch

def update_bashrc(path, content, new_version):
    new_line = f'export VYRTUOUS_VERSION="{new_version[0]}.{new_version[1]}.{new_version[2]}"'
    if "export VYRTUOUS_VERSION=" in content:
        updated = re.sub(r'export VYRTUOUS_VERSION="\d+\.\d+\.\d+"', new_line, content)
    else:
        updated = content.strip() + "\n" + new_line + "\n"
    with open(path, "w") as f:
        f.write(updated)

def update_pyproject_version(pyproject_path, new_version):
    new_version_str = f"{new_version[0]}.{new_version[1]}.{new_version[2]}"
    with open(pyproject_path, "r") as f:
        content = f.read()
    updated = re.sub(r'version\s*=\s*"\d+\.\d+\.\d+"', f'version = "{new_version_str}"', content)
    with open(pyproject_path, "w") as f:
        f.write(updated)

if __name__ == "__main__":
    bashrc_path = os.path.expanduser("~/.bashrc_vyrtuous")
    pyproject_path = os.path.join(os.getcwd(), "pyproject.toml")
    version, content = get_current_version(bashrc_path)
    new_version = increment_version(version)
    update_bashrc(bashrc_path, content, new_version)
    update_pyproject_version(pyproject_path, new_version)
    print(f"{new_version[0]}.{new_version[1]}.{new_version[2]}")
