from mshr import Rectangle, generate_mesh
from fenics import VectorFunctionSpace, FunctionSpace, Point
from fenics import DirichletBC, Expression, Constant
import numpy as np

# Parameters
mu = 1.0
dt = 0.005
rho = 1.0
n_mesh = 50


# Create mesh
domain = Rectangle(Point(0.0, 0.0), Point(1.0, 1.0))
fluid_mesh = generate_mesh(domain, n_mesh)

# Define function spaces
V = FunctionSpace(fluid_mesh, "P", 1)
bcu_wall_1 = DirichletBC(V, Constant(3), "near(x[0],0.0)")
bcu_wall_2 = DirichletBC(V, Constant(4), "near(x[0],1.0)")
bcu_wall_3 = DirichletBC(V, Constant(1), "near(x[1],0.0)")
bcu_wall_4 = DirichletBC(V, Constant(2), "near(x[1],1.0)")
bcu = [bcu_wall_1, bcu_wall_2, bcu_wall_3, bcu_wall_4]

def extract_over_mid_x(p0, x, y, n=10):
    u_list = []
    y_list = []
    for i in range(n):
        yy = y[0] + (y[1]-y[0])*i/(n-1)
        xx = 0.5*(x[0] + x[1])
        y_list.append(yy)
        u_list.append(p0(xx,yy))
    
    return u_list, y_list 



__all__ = [
    "V",
    "fluid_mesh",
    "bcu",
    "dt",
    "rho",
    "mu",
    "extract_over_mid_x",
]
