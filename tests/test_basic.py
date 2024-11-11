import ibfenics1

assert ibfenics1.cpp.add(1, 5) == 6
assert ibfenics1.cpp.__version__ == "0.0.2"


ibfenics1.mesh.ale.test_cpp()
