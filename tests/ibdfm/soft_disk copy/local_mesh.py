from ibfenics1 import Interaction
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


def kinematic_energy(u):
    return assemble(0.5 * inner(u, u) * dx)


def total_energy(disp, u):
    pass


initial_disp = Expression(("x[0]", "x[1]"), degree=1)

# Define boundary conditions for fluid solver
def calculate_fluid_boundary_conditions(V, Q):
    noslip = DirichletBC(V, (0, 0), "near(x[0],1) || near(x[0],0) || near(x[1],0)")
    upflow = DirichletBC(V, (1, 0), "near(x[1],1)")
    pinpoint = DirichletBC(Q, 0, "near(x[0],0) && near(x[1],0)", "pointwise")
    bcu = [noslip, upflow]
    bcp = [pinpoint]
    return bcu, bcp


def boundary_quardrature_rule(disp, force):
    mesh = disp.function_space().mesh()
    bdry = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
    class DirichletBoundary(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary
    bdry = MeshFunction("size_t", mesh, 1)
    DirichletBoundary().mark(bdry, 1)
    import ibfenics1
    fade = ibfenics1.cpp.FacetIntegration(mesh, bdry, 1)
    fade.fun4(disp._cpp_object, force._cpp_object)
    facets_points, facets_values = fade.fun3(disp._cpp_object, force._cpp_object)
    facets_weights = [1/n_mesh_fluid/n_mesh_fluid]*(len(facets_points)//2)
    return facets_points, facets_values, facets_weights
