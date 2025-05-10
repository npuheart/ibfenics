import os
import swanlab
import time
from loguru import logger
from fenics import *
from local_mesh import *
from ibfenics1.nssolver import IPCSSolver
from ibfenics1.nssolver import ChorinSolver
from ibfenics1.io import (
    unique_filename,
    create_xdmf_file,
    TimeManager,
    write_paramters,
    write_excel,
    write_excel_sheets,
)

nssolver = "projection"

swanlab.login(api_key="VBxEp1UBe2606KHDM9264", save=True)
swanlab.init(
    project=os.path.splitext(os.path.basename(__file__))[0],
    experiment_name=f"dt_{dt}_Ne_{Ne}_{nssolver}",
    description="二维理想瓣膜",
)

u0 = Function(Vf, name="velocity")
p0 = Function(Qf, name="pressure")
f0 = Function(Vf, name="force")

time_manager = TimeManager(T, num_steps, 200)
bcu, bcp = calculate_fluid_boundary_conditions(Vf, Qf)
# fluid_solver = IPCSSolver(u0, p0, f0, dt, nv, bcp=bcp, bcu=bcu, rho=1.0, conv=True)
fluid_solver = ChorinSolver(u0, p0, f0, dt, nv, bcp=bcp, bcu=bcu, rho=1.0, conv=True)

file_fluid_name = unique_filename(os.path.basename(__file__), "note", "/fluid.xdmf")
file_solid_name = unique_filename(os.path.basename(__file__), "note", "/solid.xdmf")
file_fluid = create_xdmf_file(fluid_mesh.mpi_comm(), file_fluid_name)
file_solid = create_xdmf_file(solid_mesh.mpi_comm(), file_solid_name)
swanlab.config['file_fluid'] = file_fluid_name
swanlab.config['file_solid'] = file_solid_name


# Solid Mechanics
# from ufl import ln
mu_s = 5.6e5
nv_s = 0.4
lambda_s = 2*mu_s*(1-nv_s)/3.0/(1-2*nv_s)
U0 = Function(Vs, name="velocity")
X0 = interpolate(Expression(("x[0]", "x[1]"), degree=2), Vs)
F0 = Function(Vs, name="force")

dVs = TestFunction(Vs)
FF = grad(X0)
L_hat = - inner(mu_s*(FF-inv(FF).T) + lambda_s*ln(det(FF))*inv(FF).T, grad(dVs))*dx
# L_hat = - inner(mu_s*(FF-inv(FF).T), grad(dVs))*dx
# 
t = 0
start_time = time.time()
for n in range(num_steps):
    t += dt
    flow_velocity.t = t
    un, pn = fluid_solver.solve(bcu, bcp)
    fluid_solver.update(un, pn)
    # Update the force
    inter.fluid_to_solid(u0._cpp_object, U0._cpp_object)
    X0.vector()[:] = X0.vector()[:] + dt * U0.vector()[:]
    inter.evaluate_current_points(X0._cpp_object)
    b = assemble(L_hat)
    F0.vector()[:] = b[:]
    inter.solid_to_fluid(f0._cpp_object, F0._cpp_object)
    if time_manager.should_output(n):
        logger.info(f"t = {t}")
        # Write to xdmf
        file_fluid.write(u0, t)
        file_fluid.write(p0, t)
        file_solid.write(X0, t)
        # SwanLab log
        swanlab.log(
            {
                "timecost": time.time() - start_time,
                "inflow":flow_velocity(0.0,0.805)[0],
                "x_displacement": X0(2.0106,0.91+1e-4)[0] - 2.0106,
                "y_displacement": X0(2.0106,0.91+1e-4)[1] - 0.91,
                "time": t, 
                "u_max": un.vector().max(), 
                "u_norm_l2":  un.vector().norm("l2")
            }, step = int(1+t/dt_minimum)
        )




















