class Vyrtuous(object):

    def __init__(self):
        self.__init__()


def test_vyrtuous():
    vyrtuous = Vyrtuous()
    assert "__init__" in dir(vyrtuous)
