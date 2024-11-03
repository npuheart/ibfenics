import os
from loguru import logger
from fenics import *
from local_mesh import *
from ibfenics.nssolver import SAVTaylorHoodSolverBDF2
from ibfenics.io import unique_filename, create_xdmf_file, TimeManager, write_paramters, write_excel


TaylorHoodSolverBDF2_1 = SAVTaylorHoodSolverBDF2.TaylorHoodSolverBDF2_1
TaylorHoodSolverBDF2_2 = SAVTaylorHoodSolverBDF2.TaylorHoodSolverBDF2_2
modified_energy        = SAVTaylorHoodSolverBDF2.modified_energy
# calculate_SAV_2        = SAVTaylorHoodSolverBDF2.calculate_SAV_2
calculate_SAV          = SAVTaylorHoodSolverBDF2.calculate_SAV

T = 5.0            
num_steps = 5000
dt = T / num_steps  
mu = nv     
rho = 1.0         

u0 = Function(V, name="velocity")
p0 = Function(Q, name="pressure")
f0 = Function(V, name="force")

time_manager = TimeManager(T, num_steps, 20)
bcu, bcp = calculate_fluid_boundary_conditions(V, Q)
navier_stokes_solver_1 = TaylorHoodSolverBDF2_1(u0, u0, p0, dt, nv)
navier_stokes_solver_2 = TaylorHoodSolverBDF2_2(u0, u0, p0, f0, dt, nv)
W = navier_stokes_solver_1.W
file_fluid_name = unique_filename(
    os.path.basename(__file__), "note", "/fluid.xdmf")
file_fluid = create_xdmf_file(fluid_mesh.mpi_comm(), file_fluid_name)

t = 0
for n in range(num_steps):
    
    t += dt
    un, pn = fluid_solver.solve(bcu, bcp)
    fluid_solver.update(un, pn)
    logger.info('u max:{0}', un.vector().max())
    if time_manager.should_output(n):
        logger.info(f"t = {t}")
        file_fluid.write(u0, t)
        file_fluid.write(p0, t)

