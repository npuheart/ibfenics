import os
from loguru import logger
from fenics import *
from local_mesh import *
from ibfenics1.nssolver import TaylorHoodSolver
from ibfenics1.io import (
    unique_filename,
    create_xdmf_file,
    TimeManager,
    write_paramters,
    write_excel,
)

u0 = Function(V, name="velocity")
p0 = Function(Q, name="pressure")
f0 = Function(V, name="force")

time_manager = TimeManager(T, num_steps, 20)
V, Q = TaylorHoodSolver.construct_function_space_bc(u0, p0)
bcu, bcp = calculate_fluid_boundary_conditions(V, Q)
fluid_solver = TaylorHoodSolver(u0, p0, f0, dt, nv)
file_fluid_name = unique_filename(os.path.basename(__file__), "note", "/fluid.xdmf")
file_fluid = create_xdmf_file(fluid_mesh.mpi_comm(), file_fluid_name)

t = 0
for n in range(num_steps):

    t += dt
    un, pn = fluid_solver.solve(bcu, bcp)
    fluid_solver.update(un, pn)
    logger.info("u max:    {0}", un.vector().max())
    logger.info("u l2 norm:{0}", un.vector().norm("l2"))
    if time_manager.should_output(n):
        logger.info(f"t = {t}")
        file_fluid.write(u0, t)
        file_fluid.write(p0, t)
