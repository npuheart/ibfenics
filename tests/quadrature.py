import ibfenics1
from fenics import *
from mshr import *

class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary


domain = Circle(Point(0.5, 0.5), 0.25)
mesh = generate_mesh(domain, 10)
bdry = MeshFunction("size_t", mesh, 1)
DirichletBoundary().mark(bdry, 1)
fade = ibfenics1.cpp.FacetIntegration(mesh, bdry, 1)
V = VectorFunctionSpace(mesh, "CG", 1)
disp = Function(V)
force = Function(V)

disp.interpolate(Expression(("x[0]","x[1]"),degree=1))
force.interpolate(Expression(("x[0]*x[0]","x[1]*x[1]"),degree=1))

fade.fun4(disp._cpp_object, force._cpp_object)
facets_points, facets_values = fade.fun3(disp._cpp_object, force._cpp_object)
facets_weights = [1]*(len(facets_points)//2)


from fenics import *
from ibfenics1.cpp import IBMesh, IBInterpolation

# 创建一个流体网格
seperations = [32,32]
points = [Point(0, 0), Point(1, 1)]
ib_mesh = IBMesh(points, seperations, 1)
fluid_mesh = ib_mesh.mesh()

# 创建一个固体网格
solid_mesh = UnitSquareMesh(10, 10)

# 创建一个流体函数空间
Vf_1 = VectorFunctionSpace(fluid_mesh, "P", 1)
uf = interpolate(Expression(("x[0]", "x[1]"), degree=2), Vf_1)
ib_mesh.build_map(uf._cpp_object)
ib_interpolation = IBInterpolation(ib_mesh, solid_mesh)

uf = Function(Vf_1)
# void points_to_fluid(Function &fluid, const std::vector<double> &values,const std::vector<double> &points, const std::vector<double> &weights)
ib_interpolation.points_to_fluid(uf._cpp_object,facets_values,facets_points,facets_weights)
File("uf.pvd") << uf

