from fenics import *
from ElasticDisk import ElasticDisk
from local_mesh import *


disp = Function(Vs, name="displacement")
velocity = Function(Vs, name="velocity")

force = Function(Vs, name="force")
disp.interpolate(initial_disp)
disp0 = Function(Vs, name="displacement0")


solid_solver = ElasticDisk(disp.function_space(), dt, nu_s, rho)

solid_solver.P_s_N
F = solid_solver.pently_force(disp)
solid_solver.solve()

File("F.pvd") << F