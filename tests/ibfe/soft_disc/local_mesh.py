from ibfenics1 import Interaction
from dolfin import *
# from mshr import *
import numpy as np

from local_parameters import *


def advance_disp_be(disp, velocity, dt):
    disp.vector()[:] = velocity.vector()[:] * dt + disp.vector()[:]


def advance_disp_bdf2(disp, disp_, velocity, dt):
    tmp = Function(disp.function_space())
    tmp.vector()[:] = disp.vector()[:]
    disp.vector()[:] = (2.0*velocity.vector()[:] * dt +
                        4.0*disp.vector()[:] - disp_.vector()[:])/3.0
    disp_.vector()[:] = tmp.vector()[:]


def calculate_volume(X):
    volume_J = assemble(det(grad(X)) * dx)
    print("体积：", volume_J)
    return volume_J


def calculate_solid_mesh(n_mesh_solid):
    solid_mesh = generate_mesh(Circle(Point(0.6, 0.5), 0.2), n_mesh_solid)
    return solid_mesh


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
N = FacetNormal(fluid_mesh)

initial_disp = Expression(("x[0]", "x[1]"), degree=1)

def calculate_fluid_boundary_conditions(V, Q):
    noslip = DirichletBC(
        V, (0, 0), "near(x[0],1) || near(x[0],0) || near(x[1],0)")
    upflow = DirichletBC(V, (1, 0), "near(x[1],1)")
    pinpoint = DirichletBC(Q, 0, "near(x[0],0) && near(x[1],0)", "pointwise")
    bcu = [noslip, upflow]
    bcp = [pinpoint]
    return bcu, bcp


# P = nu_s * (F - inv(F).T)
def total_energy(u, disp=None):
    def kinematic_energy(u):
        return assemble(0.5 * inner(u, u) * dx)

    def potential_energy(disp):
        F = grad(disp)
        tr_C = tr(F.T * F)
        J = det(F)
        return assemble(0.5*nu_s*(tr_C - 2) * dx)
    if disp is None:
        return kinematic_energy(u)
    else:
        return kinematic_energy(u) + potential_energy(disp)
