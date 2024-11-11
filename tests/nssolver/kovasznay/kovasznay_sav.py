import os
from loguru import logger
from fenics import *
from local_mesh import *
from ibfenics1.nssolver import SAVTaylorHoodSolver
from ibfenics1.io import (
    unique_filename,
    create_xdmf_file,
    TimeManager,
    write_paramters,
    write_excel,
)

stab = False
alpha = 0.1
conv = True


TaylorHoodSolver_1 = SAVTaylorHoodSolver.TaylorHoodSolver_1
TaylorHoodSolver_2 = SAVTaylorHoodSolver.TaylorHoodSolver_2
modified_energy = SAVTaylorHoodSolver.modified_energy
CAL_SAV = SAVTaylorHoodSolver.CAL_SAV_2
construct_function_space_bc = SAVTaylorHoodSolver.construct_function_space_bc

u0 = Function(V, name="velocity")
p0 = Function(Q, name="pressure")
f0 = Function(V, name="force")

time_manager = TimeManager(T, num_steps, 20)

V, Q = construct_function_space_bc(u0, p0)
bcus_1, bcps_1 = calculate_fluid_boundary_conditions(V, Q)
bcus_2, bcps_2 = calculate_fluid_boundary_conditions_sav(V, Q)
navier_stokes_solver_1 = TaylorHoodSolver_1(u0, p0, dt, nv, stab=stab, alpha=alpha)
navier_stokes_solver_2 = TaylorHoodSolver_2(
    u0, p0, f0, dt, nv, stab=stab, alpha=alpha, conv=conv
)

file_fluid_name = unique_filename(os.path.basename(__file__), "note", "/fluid.xdmf")
file_fluid = create_xdmf_file(fluid_mesh.mpi_comm(), file_fluid_name)

t = 0
for n in range(num_steps):
    t += dt
    navier_stokes_solver_1.update(u0, p0)
    navier_stokes_solver_2.update(u0, p0)
    u1, p1 = navier_stokes_solver_1.solve(bcus_1, bcps_1)
    u2, p2 = navier_stokes_solver_2.solve(bcus_2, bcps_2)
    # SAV = CAL_SAV(En, delta, dt, alpha, h, nu, u0, u1, u2, qn, N, w, p1, p2, rho)
    SAV = 1.0
    u0.vector()[:] = u1.vector()[:] + SAV * u2.vector()[:]
    p0.vector()[:] = p1.vector()[:] - SAV * p2.vector()[:]

    logger.info("u max:    {0}", u0.vector().max())
    logger.info("u l2 norm:{0}", u0.vector().norm("l2"))
    if time_manager.should_output(n):
        logger.info(f"t = {t}")
        file_fluid.write(u0, t)
        file_fluid.write(p0, t)
