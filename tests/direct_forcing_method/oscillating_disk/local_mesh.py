from ibfenics import Interaction
from dolfin import *
from mshr import *
import numpy as np

from local_parameters import *

#更新前一步的固体位移
def advance_disp_be(disp, velocity, dt):
    disp.vector()[:] = velocity.vector()[:] * dt + disp.vector()[:]

#计算固体的体积变化
def calculate_volume(X):
    volume_J = assemble(det(grad(X)) * dx)
    print("体积：", volume_J)
    return volume_J

#定义固体网格
def calculate_solid_mesh(n_mesh_solid):
    solid_mesh = generate_mesh(Circle(Point(0.5, 0.5), 0.2), n_mesh_solid)
    return solid_mesh

#创建固体网格
solid_mesh = calculate_solid_mesh(n_mesh_solid)

orders = [order_velocity, order_pressure, order_displacement]
seperations = [n_mesh_fluid, n_mesh_fluid]
box_points = [Point(0, 0), Point(1, 1)]
interaction = Interaction(box_points, seperations, solid_mesh, orders)

fluid_mesh = interaction.fluid_mesh
ib_mesh = interaction.ib_mesh
ib_interpolation = interaction.ib_interpolation
Vs = interaction.Vs
Vf = interaction.Vf
Vf_1 = interaction.Vf_1
Vp = interaction.Vp


def kinematic_energy(u):
    return assemble(0.5 * inner(u, u) * dx)


def total_energy(disp, u):
    pass

initial_disp = Expression(("x[0]", "x[1]"), degree=1)
initial_velocity = Expression(('psi0*2*pi*sin(2*pi*x[0])*cos(2*pi*x[1])', '-psi0*2*pi*cos(2*pi*x[0])*sin(2*pi*x[1])'), degree=2, psi0=0.05)
boundary_velocity = Expression(('psi0*2*pi*sin(2*pi*x[0])*cos(2*pi*x[1])', '-psi0*2*pi*cos(2*pi*x[0])*sin(2*pi*x[1])'), degree=2, psi0=0.05)
# Define boundary conditions for fluid solver
def calculate_fluid_boundary_conditions(V, Q):
    noslip = DirichletBC(V, boundary_velocity, "near(x[0],1) || near(x[0],0) || near(x[1],0)|| near(x[1],1)")
    # pinpoint = DirichletBC(Q, 0, "near(x[0],0) && near(x[1],0)", "pointwise")
    bcu = [noslip]
    bcp = []
    return bcu, bcp

