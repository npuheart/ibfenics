from ibfe import cpp


def test_cpp():
    assert cpp.add(1,5) == 7
    assert cpp.__version__ == "0.0.2"



