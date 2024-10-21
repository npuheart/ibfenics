



from fenics import *
from mshr import *

def get_mesh(n_mesh_solid = 40):
    circle_outer = Circle(Point(0.5,0.5), 0.25)
    circle_inner = Circle(Point(0.5,0.5), 0.25-0.0625)
    solid_mesh   = generate_mesh(circle_outer-circle_inner, n_mesh_solid)
    return solid_mesh
