def test_read():
    read_the_docs()


def read_the_docs():
    with open("tests/unit_tests/docs.python.org", "r") as f:
        print(f.read())


if __name__ == "__main__":
    test_read()
