from fenics import *
from mshr import *
from ibfenics1 import Interaction


def calculate_solid_mesh(n_mesh_solid=40):
    circle_outer = Circle(Point(0.5, 0.5), 0.25 + 0.0625)
    circle_inner = Circle(Point(0.5, 0.5), 0.25)
    solid_mesh = generate_mesh(circle_outer - circle_inner, n_mesh_solid)
    return solid_mesh


# 32 40
# 64 80
# 128 160
n_mesh_fluid = 32
n_mesh_solid = 40

solid_mesh = calculate_solid_mesh(n_mesh_solid)
order_velocity = 2
order_pressure = 1
order_displacement = 1

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

# 暴露标识符
__all__ = ["solid_mesh", "fluid_mesh", "Vs", "Vf", "Vf_1", "Vp", "ib_interpolation"]
