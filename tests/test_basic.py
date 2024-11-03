import ibfenics

assert ibfenics.cpp.add(1, 5) == 6
assert ibfenics.cpp.__version__ == "0.0.2"


ibfenics.mesh.ale.test_cpp()
