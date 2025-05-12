import os
import swanlab
import time
from loguru import logger
from fenics import *
from ufl import cofac
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
CC = 100000.0  # dyn/cm^2 (2kPa)
bf = 1.0
bt = 1.0
bfs = 1.0
kappa = 1e6 # Penalty for incompressibility
beta  = 5e7  # Penalty for boundary conditions

pressure = Expression("t", degree=2, t=0.0)
ds = Measure('ds', domain=solid_mesh, subdomain_data=boundaries)

nssolver = "projection"

swanlab.login(api_key="VBxEp1UBe2606KHDM9264", save=True)
swanlab.init(
    project=os.path.splitext(os.path.basename(__file__))[0],
    experiment_name=f"dt_{dt}_Ne_{Ne}_Nl_{Nl}",
    description="三维理想左心室",
    config={'dt': dt, 'Ne': Ne, 'Nl': Nl},
)

u0 = Function(Vf, name="velocity")
p0 = Function(Qf, name="pressure")
f0 = Function(Vf, name="force")

time_manager = TimeManager(T, num_steps, 200)
bcu, bcp = calculate_fluid_boundary_conditions(Vf, Qf)
fluid_solver = ChorinSolver(u0, p0, f0, dt, nv, bcp=bcp, bcu=bcu, rho=1.0, conv=True)

file_fluid_name = unique_filename(os.path.basename(__file__), "note", "/fluid.xdmf")
file_solid_name = unique_filename(os.path.basename(__file__), "note", "/solid.xdmf")
file_fluid = create_xdmf_file(fluid_mesh.mpi_comm(), file_fluid_name)
file_solid = create_xdmf_file(solid_mesh.mpi_comm(), file_solid_name)
swanlab.config['file_fluid'] = file_fluid_name
swanlab.config['file_solid'] = file_solid_name


# Solid Mechanics
# from ufl import ln
mu_s = 0.1
nv_s = 0.4
lambda_s = 2*mu_s*(1-nv_s)/3.0/(1-2*nv_s)
U0 = Function(Vs, name="velocity")
X0 = interpolate(Expression(("x[0]", "x[1]", "x[2]"), degree=2), Vs)
X_start = interpolate(Expression(("x[0]", "x[1]", "x[2]"), degree=2), Vs)
F0 = Function(Vs, name="force")

dVs = TestFunction(Vs)
FF = grad(X0)

def first_PK_stress(X0):                           # X is the position of current configuration
    I = Identity(len(X0))
    F = variable(grad(X0))                         # nabla_grad is used someplaces. I think grad is correct.
    J = det(F)
    C = pow(J, -float(2)/3) * F.T*F
    E = 0.5*(C - I)
    E11, E12, E13 = inner(E*e1, e1), inner(E*e1, e2), inner(E*e1, e3)
    E21, E22, E23 = inner(E*e2, e1), inner(E*e2, e2), inner(E*e2, e3)
    E31, E32, E33 = inner(E*e3, e1), inner(E*e3, e2), inner(E*e3, e3)
    Q = bf*E11**2 + bt*(E22**2 + E33**2 + E23**2 + E32**2) \
      + bfs*(E12**2 + E21**2 + E13**2 + E31**2)
    # passive strain energy
    Wpassive = CC/2.0 * (exp(Q) - 1)
    # incompressibility
    Winc = kappa * ln(J)**2
    return diff(Wpassive+Winc, F)

L_hat = - inner(first_PK_stress(X0), grad(dVs))*dx
L_hat += beta*inner(X_start-X0, TestFunction(Vs))*ds(5)
L_hat -= pressure * inner(cofac(FF)*FacetNormal(solid_mesh), TestFunction(Vs))*ds(6)

# 
t = 0
start_time = time.time()
for n in range(num_steps):
    t += dt
    pressure.t = t
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
                "time": t, 
                "u_max": un.vector().max(), 
                "u_norm_l2":  un.vector().norm("l2")
            }, step = int(1+t/dt_minimum)
        )















