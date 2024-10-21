
# 导入网格
import os
from fenics import *

marker_circle = 1
marker_beam = 2

 
mesh_path = os.path.expanduser("~") + "/mesh/benchmark/circle_beam/"

solid_mesh = Mesh()
with XDMFFile(mesh_path + "circle_beam_40.xdmf") as xdmf:
    xdmf.read(solid_mesh)

bdry = MeshFunction("size_t", solid_mesh, mesh_path + "circle_beam_40_boundaries.xml")
domains = MeshFunction("size_t", solid_mesh, mesh_path +  "circle_beam_40_domains.xml")
dx = Measure("dx")(subdomain_data=domains)

if __name__ == "__main__":
    File("bdry.pvd") << bdry
    File("domains.pvd") << domains
    File("solid_mesh.pvd") << solid_mesh