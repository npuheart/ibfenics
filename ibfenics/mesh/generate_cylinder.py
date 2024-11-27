from ibfenics.io import (
    unique_filename,
    create_xdmf_file,
    write_excel,
    write_bg_mesh,
    write_mesh,
)
from ibfenics.IU import IU
from dolfin import *
from mshr import *
import os

mesh_index = 9999
mesh_name = "cylinder"
mesh_path = os.path.expanduser("~") + "/mesh/benchmark/cylinder/"
mesh_size = 40

L = 2.20 * IU.m
H = 0.41 * IU.m
r = 0.05 * IU.m


origin = Point(0.0 * IU.m, 0.0 * IU.m)  # point O
center = Point(0.2 * IU.m, 0.2 * IU.m)  # point C


_domain_1 = Circle(center, r)       #固体圆柱
_domain = _domain_1 
_domain.set_subdomain(1, _domain_1)

# output mesh
_mesh = generate_mesh(_domain, mesh_size)

# output domains
_domains = MeshFunction("size_t", _mesh, 2, _mesh.domains())

# output boundaries


class Boundary_2(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary


boundary = MeshFunction("size_t", _mesh, 1)
Boundary_2().mark(boundary, 4)

# write mesh
write_mesh(_mesh, _domains, boundary, mesh_size, mesh_name, mesh_path)

# write background mesh
write_bg_mesh(L, H, mesh_size, mesh_name, mesh_path)
