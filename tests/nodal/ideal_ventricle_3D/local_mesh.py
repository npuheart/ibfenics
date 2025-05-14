from fenics import VectorFunctionSpace, FunctionSpace, Point
from fenics import DirichletBC, Expression, Constant, interpolate, File
from fenics import BoxMesh, Mesh, MeshEditor, CellType, MeshFunction, Function, XDMFFile
from ibfenics1.cpp import IBMesh3D, IBInterpolation3D
import numpy as np

def read_fibers(mesh, data_file, W, index = 0):
    f0 = Function(W)
    s0 = Function(W)
    n0 = Function(W)
    f_in = XDMFFile(mesh.mpi_comm(), data_file)
    f_in.read_checkpoint(f0, "fiber", 0)
    f_in.read_checkpoint(s0, "sheet", 0)
    f_in.read_checkpoint(n0, "normal", 0)
    f_in.close()
    return f0,s0,n0

def read_mesh_general(mesh_path):
    mesh_file = XDMFFile(mesh_path)
    mesh = Mesh()
    mesh_file.read(mesh)
    mesh_file.close()
    return mesh



# Parameters
mesh_path = "/kokkos_2/mesh-benchmark/B009-ldrb/"
nv = 0.04
T = 1.0
dt = 1e-5
num_steps = int(T/dt)
rho = 1.0
Nl = 4
Ne = 20
dt_minimum = 1e-5

# Mesh
order_velocity = 2
order_pressure = 1
order_displacement = 2


ib_mesh = IBMesh3D([Point(0, 0, 0), Point(5, 5, 5)], [Ne, Ne, Ne], order_velocity)
fluid_mesh = ib_mesh.mesh()
solid_mesh = read_mesh_general(mesh_path+"mesh.xdmf")
boundaries = MeshFunction("size_t", solid_mesh, mesh_path + "boundaries.xml")


for i in range(len(solid_mesh.coordinates())):
    solid_mesh.coordinates()[i][0] *= 0.1
    solid_mesh.coordinates()[i][1] *= 0.1
    solid_mesh.coordinates()[i][2] *= 0.1
for i in range(len(solid_mesh.coordinates())):
    solid_mesh.coordinates()[i][0] += 3.0
    solid_mesh.coordinates()[i][1] += 2.5
    solid_mesh.coordinates()[i][2] += 2.5
        
# Define function spaces
Qf = FunctionSpace(fluid_mesh, "P", order_pressure)
Vf = VectorFunctionSpace(fluid_mesh, "P", order_velocity)
Vs = VectorFunctionSpace(solid_mesh, "P", order_displacement)

e1, e2, e3 = read_fibers(solid_mesh, mesh_path+"microstructure.xdmf", Vs)


# 初始化
uf = interpolate(Expression(("x[0]", "x[1]", "x[2]"), degree=2), Vf)
ib_mesh.build_map(uf._cpp_object)
inter = IBInterpolation3D(ib_mesh, solid_mesh)
us = interpolate(Expression(("x[0]", "x[1]", "x[2]"), degree=2), Vs)
inter.evaluate_current_points(us._cpp_object)



File("bdry.pvd") << boundaries
File("mesh.pvd") << solid_mesh
# File("fibers.pvd") << f0
# File("sheet.pvd") << s0
# File("normal.pvd") << n0

def calculate_fluid_boundary_conditions(Vf, Qf):
    flow_velocity = Expression(("0.0", "0.0", "0.0"), degree=1)
    bcu_inflow = DirichletBC(Vf, flow_velocity, "near(x[1],5.0)")
    bcu_wall = DirichletBC(Vf, Expression(("0.0", "0.0", "0.0"), degree=1), 
                           "near(x[0],0.0) || near(x[0],5.0) || near(x[1],0.0) || near(x[2],0.0) || near(x[2],5.0)")
    # bcp_outlet = DirichletBC(Qf, Expression("0.0", degree=1), "near(x[0],8.0)")
    bcu = [bcu_inflow, bcu_wall]
    bcp = []
    return bcu, bcp


def extract_over_mid_y(u0, p0, x, y, n=10):
    u_list = []
    v_list = []
    p_list = []
    x_list = []
    for i in range(n):
        xx = x[0] + (x[1]-x[0])*i/(n-1)
        yy = 0.5*(y[0] + y[1])
        x_list.append(xx)
        p_list.append(p0(xx,yy))
        u_list.append(u0(xx,yy)[0])
        v_list.append(u0(xx,yy)[1])
    
    return u_list, v_list, p_list, x_list 

def extract_over_mid_x(u0, p0, x, y, n=10):
    u_list = []
    v_list = []
    p_list = []
    y_list = []
    for i in range(n):
        xx = 0.5*(x[0] + x[1])
        yy = y[0] + (y[1]-y[0])*i/(n-1)
        y_list.append(yy)
        p_list.append(p0(xx,yy))
        u_list.append(u0(xx,yy)[0])
        v_list.append(u0(xx,yy)[1])
    
    return u_list, v_list, p_list, y_list 


__all__ = [
    "Vs",
    "Vf",
    "Qf",
    "e1","e2","e3",
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
    "extract_over_mid_x",
    "extract_over_mid_y",
]
