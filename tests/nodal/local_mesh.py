from fenics import VectorFunctionSpace, FunctionSpace, Point
from fenics import DirichletBC, Expression, Constant
from fenics import RectangleMesh
import numpy as np

# Parameters
nv = 0.001
T = 1.0
dt = 0.005
num_steps = int(T/dt)
rho = 1.0
n_mesh = 50
dt_minimum = 1e-4

# Create mesh
fluid_mesh = RectangleMesh(Point(0.0, 0.0), Point(1.0, 1.0), n_mesh, n_mesh)

# Define function spaces
V = VectorFunctionSpace(fluid_mesh, "P", 2)
Q = FunctionSpace(fluid_mesh, "P", 1)

def calculate_fluid_boundary_conditions(V, Q):
    flow_velocity = Expression(("1.0", "0.0"), degree=1)
    bcu_inflow = DirichletBC(V, flow_velocity, "near(x[1],1.0)")
    bcu_wall = DirichletBC(V, Expression(("0.0", "0.0"), degree=1), "near(x[1],0.0) || near(x[0],0.0) || near(x[0],1.0)")
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
    "V",
    "Q",
    "fluid_mesh",
    "calculate_fluid_boundary_conditions",
    "nv",
    "T",
    "num_steps",
    "dt_minimum",
    "n_mesh",
    "dt",
    "rho",
    "extract_over_mid_x",
    "extract_over_mid_y",
]
