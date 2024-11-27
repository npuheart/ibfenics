import os
from ibfenics import Interaction
from dolfin import *
from mshr import *
import numpy as np

from local_parameters import *


def advance_disp_be(disp, velocity, dt):
    disp.vector()[:] = velocity.vector()[:] * dt + disp.vector()[:]


def calculate_volume(X):
    volume_J = assemble(det(grad(X)) * dx)
    print("体积：", volume_J)
    return volume_J


marker_circle = 1
mesh_path = os.path.expanduser("~") + "/mesh/benchmark/circle_beam/"
def calculate_solid_mesh(n_mesh_solid):
    solid_mesh = Mesh()
    with XDMFFile(mesh_path + "circle_beam_40.xdmf") as xdmf:
        xdmf.read(solid_mesh)
    return solid_mesh


solid_mesh = calculate_solid_mesh(n_mesh_solid)
bdry = MeshFunction("size_t", solid_mesh, mesh_path + "circle_beam_40_boundaries.xml")
domains = MeshFunction("size_t", solid_mesh, mesh_path + "circle_beam_40_domains.xml")
ddx = Measure("dx")(subdomain_data=domains)


orders = [order_velocity, order_pressure, order_displacement]
seperations = [n_mesh_fluid * 6, n_mesh_fluid]
box_points = [Point(0, 0), Point(2.2, 0.41)]
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

# Define boundary conditions for fluid solver
def calculate_fluid_boundary_conditions(V, Q):
    inflow = "near(x[0], 0)"
    outflow = "near(x[0], 2.2)"
    walls = "near(x[1], 0) || near(x[1], 0.41)"
    cylinder = "on_boundary && x[0]>0.1 && x[0]<0.3 && x[1]>0.1 && x[1]<0.3"
    inflow_profile = ("4.0*1.5*x[1]*(0.41 - x[1]) / pow(0.41, 2)", "0")
    bcu_inflow = DirichletBC(V, Expression(inflow_profile, degree=2), inflow)
    bcu_walls = DirichletBC(V, Constant((0, 0)), walls)
    bcu_cylinder = DirichletBC(V, Constant((0, 0)), cylinder)
    bcp_outflow = DirichletBC(Q, Constant(0), outflow)
    bcus = [bcu_inflow, bcu_walls, bcu_cylinder]
    bcps = [bcp_outflow]
    return bcus, bcps



if __name__ == "__main__":
    File("bdry.pvd") << bdry
    File("domains.pvd") << domains
    File("solid_mesh.pvd") << solid_mesh
