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
    circle_outer = Circle(Point(0.5, 0.5), 0.25 + 0.0625)
    circle_inner = Circle(Point(0.5, 0.5), 0.25)
    solid_mesh = generate_mesh(circle_outer - circle_inner, n_mesh_solid)
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


def kinematic_energy(u):
    return assemble(0.5 * inner(u, u) * dx)


def total_energy(disp, u):
    pass


class InitialDisplacement(UserExpression):
    def eval(self, value, x):
        R = 0.25
        gamma = 0.15
        center = [0.5, 0.5]
        s = [0.0, 0.0]
        s[1] = (
            np.sqrt(
                (x[1] - center[1]) * (x[1] - center[1])
                + (x[0] - center[0]) * (x[0] - center[0])
            )
            - R
        )
        s[0] = R * np.arccos((x[0] - 0.5) / (R + s[1]))

        if x[1] < 0.5:
            s[0] = 2.0 * np.pi * R - s[0]
        else:
            s[0] = s[0]
        value[0] = (R + s[1]) * np.cos(s[0] / R) + 0.5
        value[1] = (R + s[1] + gamma) * np.sin(s[0] / R) + 0.5

    def value_shape(self):
        return (2,)


initial_disp = InitialDisplacement()

def calculate_fluid_boundary_conditions(V, Q):
    noslip = DirichletBC(
        V, (0, 0), "near(x[0],1) || near(x[0],0) || near(x[1],0)")
    upflow = DirichletBC(V, (0, 0), "near(x[1],1)")
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
        return assemble(0.5*nu_s*(tr_C - J) * dx)
        # return assemble(0.5*nu_s*(tr_C - 2.0 - ln(J)) * dx)
    if disp is None:
        return kinematic_energy(u)
    else:
        return kinematic_energy(u) + potential_energy(disp)
