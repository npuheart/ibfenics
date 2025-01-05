import ibfenics1
from fenics import *

class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary


mesh = UnitSquareMesh(10,10)
bdry = MeshFunction("size_t", mesh, 1)
DirichletBoundary().mark(bdry, 1)
fade = ibfenics1.cpp.FacetIntegration(mesh, bdry,1)
V = FunctionSpace(mesh, "CG", 1)
disp = Function(V)
force = Function(V)
fade.fun4(disp._cpp_object, force._cpp_object)


