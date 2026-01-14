import os
import re
import sys


def find_channel_snowflake_in_triple_sql(root):
    extensions = {".py"}
    results = []
    sql_block_pattern = re.compile(r"'''([\s\S]*?)'''", re.MULTILINE)
    target_pattern = re.compile(r"\bchannel_snowflake\b", re.IGNORECASE)
    for dirpath, _, filenames in os.walk(root):
        for filename in filenames:
            if os.path.splitext(filename)[1] in extensions:
                path = os.path.join(dirpath, filename)
                try:
                    with open(path, "r", encoding="utf-8") as f:
                        content = f.read()
                    for block in sql_block_pattern.finditer(content):
                        sql_text = block.group(1)
                        if target_pattern.search(sql_text):
                            line_number = content[: block.start()].count("\n") + 1
                            results.append((path, line_number))
                except Exception:
                    pass
    return results


def main():
    root = sys.argv[1] if len(sys.argv) > 1 else "."
    results = find_channel_snowflake_in_triple_sql(root)
    for path, line in results:
        print(f"{path}:{line}")
    print(f"\nTotal matches: {len(results)}")


if __name__ == "__main__":
    main()
