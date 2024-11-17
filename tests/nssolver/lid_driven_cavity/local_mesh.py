from mshr import Rectangle, generate_mesh
from fenics import VectorFunctionSpace, FunctionSpace, Point
from fenics import DirichletBC, Expression, Constant
import numpy as np

# Parameters
nv = 1.0
T = 0.5
num_steps = 100
dt = T / num_steps
rho = 1.0
n_mesh = 50


# Create mesh
domain = Rectangle(Point(0.0, 0.0), Point(1.0, 1.0))
fluid_mesh = generate_mesh(domain, n_mesh)

# Define function spaces
V = VectorFunctionSpace(fluid_mesh, "P", 2)
Q = FunctionSpace(fluid_mesh, "P", 1)


def calculate_fluid_boundary_conditions(V, Q):
    # Define boundaries
    all_boundary = "on_boundary"

    flow_velocity = Expression(("1.0", "0.0"), degree=1)

    # Define boundary conditions
    bcu_inflow = DirichletBC(V, flow_velocity, "near(x[1],1.0)")
    bcu_wall = DirichletBC(V, Expression(("0.0", "0.0"), degree=1), "near(x[1],0.0) || near(x[0],0.0) || near(x[0],1.0)")
    bcp_1 = DirichletBC(Q, Constant(0), "near(x[1],0.0) && near(x[0],0.0)", "pointwise")
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
    "dt",
    "rho",
    "extract_over_mid_x",
    "extract_over_mid_y",
]
