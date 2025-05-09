from fenics import VectorFunctionSpace, FunctionSpace, Point
from fenics import DirichletBC, Expression, Constant, interpolate
from fenics import RectangleMesh
from ibfenics1.cpp import IBMesh, IBInterpolation
import numpy as np

# Parameters
nv = 0.1
T = 10.0
dt = 0.00005
num_steps = int(T/dt)
rho = 1.0
Nl = 4
Ne = 32
dt_minimum = 1e-5

# Mesh
order_velocity = 2
order_pressure = 1
order_displacement = 2

ib_mesh = IBMesh([Point(0, 0), Point(8, 1.61)], [4*Ne, Ne], order_velocity)
fluid_mesh = ib_mesh.mesh()
solid_mesh = RectangleMesh(Point(2.0, 0.0), Point(2.0212, 0.7), Nl, Nl*10)

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











flow_velocity = Expression(("5 * (sin(2 * M_PI * t) + 1.1) * x[1] * (1.61 - x[1]);", "0.0"), degree=1, t=0.0)


def calculate_fluid_boundary_conditions(Vf, Qf):
    bcu_inflow = DirichletBC(Vf, flow_velocity, "near(x[0],0.0)")
    bcu_wall = DirichletBC(Vf, Expression(("0.0", "0.0"), degree=1), "near(x[1],0.0) || near(x[1],1.61)")
    bcp_outlet = DirichletBC(Qf, Expression("0.0", degree=1), "near(x[0],8.0)")
    bcu = [bcu_inflow, bcu_wall]
    bcp = [bcp_outlet]
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
    "flow_velocity",
    "inter",
    "fluid_mesh",
    "solid_mesh",
    "calculate_fluid_boundary_conditions",
    "nv",
    "T",
    "num_steps",
    "dt_minimum",
    "Ne",
    "dt",
    "rho",
    "extract_over_mid_x",
    "extract_over_mid_y",
]
