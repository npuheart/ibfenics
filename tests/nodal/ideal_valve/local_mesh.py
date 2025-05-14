from fenics import VectorFunctionSpace, FunctionSpace, Point
from fenics import DirichletBC, Expression, Constant, interpolate
from fenics import RectangleMesh, Mesh, MeshEditor, CellType, near, MeshFunction, SubDomain
from ibfenics1.cpp import IBMesh, IBInterpolation
import numpy as np

def connect_mesh(mesh0, mesh1):
    mesh = Mesh()
    editor = MeshEditor()
    mesh0_vertices = mesh0.num_vertices()
    mesh1_vertices = mesh1.num_vertices()
    mesh0_cells = mesh0.num_cells()
    mesh1_cells = mesh1.num_cells()
    editor.open(mesh, "triangle", 2, 2)
    editor.init_vertices(mesh0_vertices + mesh1_vertices)
    editor.init_cells(mesh0_cells + mesh1_cells)
    # 先把第一个网格的单元移动到新的网格上
    for i in range(mesh0_vertices):
        point = mesh0.coordinates()[i]
        node = Point(point[0], point[1])
        editor.add_vertex(i, node)
    for i in range(mesh1_vertices):
        point = mesh1.coordinates()[i]
        node = Point(point[0], point[1])
        editor.add_vertex(i+mesh0_vertices, node)
    # 把第一个网格的单元移动到新的网格上
    for i in range(mesh0_cells):
        cell = mesh0.cells()[i]
        element = [cell[0], cell[1], cell[2]]
        editor.add_cell(i, element)
    for i in range(mesh1_cells):
        cell = mesh1.cells()[i]
        element = [cell[0]+mesh0_vertices, cell[1]+mesh0_vertices, cell[2]+mesh0_vertices]
        editor.add_cell(i+mesh0_cells, element)
    editor.close()
    return mesh

# Parameters
nv = 10.0
T = 3.0
dt = 5e-5
num_steps = int(T/dt)
rho = 100.0
Nl = 4
Ne = 40
dt_minimum = 1e-5

# Mesh
order_velocity = 2
order_pressure = 1
order_displacement = 2

ib_mesh = IBMesh([Point(0, 0), Point(8, 1.61)], [4*Ne, Ne], order_velocity)
fluid_mesh = ib_mesh.mesh()
solid_mesh_0 = RectangleMesh(Point(2.0-0.0212, 0.0), Point(2.0, 0.7), Nl, Nl*10)
# solid_mesh = RectangleMesh.create([Point(2.0, 0.0), Point(2.0212, 0.7)], [Nl, Nl*20], CellType.Type.quadrilateral)
solid_mesh_1 = RectangleMesh(Point(2.0-0.0212, 0.91), Point(2.0, 1.61), Nl, Nl*10)
solid_mesh = connect_mesh(solid_mesh_0, solid_mesh_1)

boundaries = MeshFunction("size_t", solid_mesh, solid_mesh.topology().dim()-1)
boundaries.set_all(0)
class BottomBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return near(x[1], 0.0) or near(x[1], 1.61) 

bottom = BottomBoundary()
bottom.mark(boundaries, 1)  

# Define function spaces
Qf = FunctionSpace(fluid_mesh, "P", order_pressure)
Vf = VectorFunctionSpace(fluid_mesh, "P", order_velocity)
Vs = VectorFunctionSpace(solid_mesh, "P", order_displacement)

# 初始化
uf = interpolate(Expression(("x[0]", "x[1]"), degree=2), Vf)
ib_mesh.build_map(uf._cpp_object)
inter = IBInterpolation(ib_mesh, solid_mesh)
us = interpolate(Expression(("x[0]", "x[1]"), degree=2), Vs)
inter.evaluate_current_points(us._cpp_object)










flow_velocity = Expression(("5 * (sin(2 * M_PI * t) + 1.1) * x[1] * (1.61 - x[1]);", "0.0"), degree=2, t=0.0)


def calculate_fluid_boundary_conditions(Vf, Qf):
    bcu_inflow = DirichletBC(Vf, flow_velocity, "near(x[0],0.0)")
    bcu_wall = DirichletBC(Vf, Expression(("0.0", "0.0"), degree=1), "near(x[1],0.0) || near(x[1],1.61)")
    bcp_outlet = DirichletBC(Qf, Expression("0.0", degree=1), "near(x[0],8.0)")
    bcu = [bcu_inflow, bcu_wall]
    bcp = [bcp_outlet]
    return bcu, bcp


__all__ = [
    "Vs",
    "Vf",
    "Qf",
    "flow_velocity",
    "inter",
    "fluid_mesh",
    "solid_mesh",
    "boundaries",
    "calculate_fluid_boundary_conditions",
    "nv",
    "T",
    "num_steps",
    "dt_minimum",
    "Ne",
    "Nl",
    "dt",
    "rho",
]
