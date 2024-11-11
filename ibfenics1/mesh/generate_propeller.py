from ibfenics1.io import (
    unique_filename,
    create_xdmf_file,
    write_excel,
    write_bg_mesh,
    write_mesh,
)
from ibfenics1.IU import IU
from dolfin import *
from mshr import *
import os

mesh_index = 9998
mesh_name = "propeller"
mesh_path = os.path.expanduser("~") + "/mesh/benchmark/propeller/"
mesh_size = 40

L = 1.00 * IU.m
H = 1.00 * IU.m
r = 0.05 * IU.m
l = 0.25 * IU.m
h = 0.02 * IU.m

origin = Point(0.0 * IU.m, 0.0 * IU.m)  # point O
center = Point(0.5 * IU.m, 0.5 * IU.m)  # point C
beam_1_point_1 = Point(center.x() + r + l, center.y() + r)
beam_1_point_2 = Point(center.x() - r - l, center.y() - r)
beam_2_point_1 = Point(center.x() + r, center.y() + l + r)
beam_2_point_2 = Point(center.x() - r, center.y() - l - r)

_domain_1 = Rectangle(beam_1_point_1, beam_1_point_2)
_domain_2 = Rectangle(beam_2_point_1, beam_2_point_2)
_domain = _domain_1 + _domain_2

# output mesh
_mesh = generate_mesh(_domain, mesh_size)

# output domains
_domains = MeshFunction("size_t", _mesh, 2, _mesh.domains())

# output boundaries
class Boundary_1(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < center.x() + r and on_boundary


class Boundary_2(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary


boundary = MeshFunction("size_t", _mesh, 1)
Boundary_2().mark(boundary, 4)
Boundary_1().mark(boundary, 3)

# write mesh
write_mesh(_mesh, _domains, boundary, mesh_size, mesh_name, mesh_path)

# write background mesh
write_bg_mesh(L, H, mesh_size, mesh_name, mesh_path)
