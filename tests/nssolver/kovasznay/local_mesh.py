from mshr import Rectangle, generate_mesh
from fenics import VectorFunctionSpace, FunctionSpace, Point
from fenics import DirichletBC, Expression, Constant
import numpy as np

# Parameters
nv = 0.025
T = 5.0
num_steps = 5000
dt = T / num_steps
rho = 1.0


# Create mesh
domain = Rectangle(Point(-0.5, -0.5), Point(1.0, 0.5))
fluid_mesh = generate_mesh(domain, 64)

# Define function spaces
V = VectorFunctionSpace(fluid_mesh, "P", 2)
Q = FunctionSpace(fluid_mesh, "P", 1)


def calculate_fluid_boundary_conditions(V, Q):
    # Define boundaries
    all_boundary = "on_boundary"

    # Define inflow profile
    lam = 0.5 / nv - np.sqrt(0.25 / nv ** 2 + 4 * np.pi ** 2)
    flow_velocity = Expression(
        ("1-exp(lam*x[0])*cos(2*pi*x[1])", "lam/2.0/pi*exp(lam*x[0])*sin(2*pi*x[1])"),
        degree=2,
        lam=lam,
        pi=np.pi,
    )

    # Define boundary conditions
    bcu_inflow = DirichletBC(V, flow_velocity, all_boundary)
    bcp_1 = DirichletBC(Q, Constant(0), "near(x[1],0.0) && near(x[0],0.0)", "pointwise")
    bcu = [bcu_inflow]
    bcp = [bcp_1]
    return bcu, bcp


def calculate_fluid_boundary_conditions_sav(V, Q):
    bcu_3 = DirichletBC(V, Constant((0, 0)), "on_boundary")
    bcp_1 = DirichletBC(Q, Constant(0), "near(x[1],0.0) && near(x[0],0.0)", "pointwise")
    bcus_2 = [bcu_3]
    bcps_2 = [bcp_1]
    return bcus_2, bcps_2


__all__ = [
    "V",
    "Q",
    "fluid_mesh",
    "calculate_fluid_boundary_conditions",
    "calculate_fluid_boundary_conditions_sav",
    "nv",
    "T",
    "num_steps",
    "dt",
    "rho",
]
