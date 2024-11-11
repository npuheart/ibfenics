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

mesh_index = 9999
mesh_name = "circle_beam"
mesh_path = os.path.expanduser("~") + "/mesh/benchmark/circle_beam/"
mesh_size = 80

L = 2.5 * IU.m
H = 0.41 * IU.m
r = 0.05 * IU.m
l = 0.35 * IU.m
h = 0.02 * IU.m

origin = Point(0.0 * IU.m, 0.0 * IU.m)  # point O
center = Point(0.2 * IU.m, 0.2 * IU.m)  # point C
end = Point(center.x() + r + l, center.y())  # point A

beam_1_point = Point(center.x(), center.y() - 0.5 * h)
beam_2_point = Point(end.x(), end.y() + 0.5 * h)

_domain_1 = Circle(center, r)
_domain_2 = Rectangle(beam_1_point, beam_2_point)
_domain = _domain_1 + _domain_2
_domain.set_subdomain(2, _domain_2)
_domain.set_subdomain(1, _domain_1)

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
