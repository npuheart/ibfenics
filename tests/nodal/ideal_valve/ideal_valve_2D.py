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
beta = 1e8
# C0=2e5
C0=0
C1=1e6
# kappa_s=4e5
kappa_s=0

ds = Measure('ds', domain=solid_mesh, subdomain_data=boundaries)

nssolver = "projection"
material = "f1" # "n1", "n2", "f0", "f1", "f2", "f3"

swanlab.login(api_key="VBxEp1UBe2606KHDM9264", save=True)
swanlab.init(
    project=os.path.splitext(os.path.basename(__file__))[0],
    experiment_name=f"dt_{dt}_beta_{beta}_C0_{C0}_C1_{C1}_ks_{kappa_s}_Ne_{Ne}_Nl_{Nl}_{material}",
    description="二维理想瓣膜",
    config={'dt': dt, 'beta': beta, 'C0': C0, 'C1':C1, 'ks':kappa_s, 'Ne': Ne, 'Nl': Nl, 'material': material},
    
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

U0 = Function(Vs, name="velocity")
X0 = interpolate(Expression(("x[0]", "x[1]"), degree=2), Vs)
X_start = interpolate(Expression(("x[0]", "x[1]"), degree=2), Vs)
F0 = Function(Vs, name="force")

def MaterialModel_1(X0):
    # 计算形变梯度
    FF = grad(X0)
    # 本构参数
    mu_s = 5.6e5
    nv_s = 0.4
    lambda_s = 2*mu_s*(1-nv_s)/3.0/(1-2*nv_s)
    dVs = TestFunction(Vs)
    L_hat = - inner(mu_s*(FF-inv(FF).T) + lambda_s*ln(det(FF))*inv(FF).T, grad(dVs))*dx
    return L_hat

def MaterialModel_2(X0):
    # 计算形变梯度
    FF = grad(X0)
    # 本构参数
    mu_s = 5.6e5
    dVs = TestFunction(Vs)
    L_hat = - inner(mu_s*(FF-inv(FF).T), grad(dVs))*dx
    return L_hat

f1 = Constant((1.0, 0.0))
f2 = Constant((0.7071067811865475,0.7071067811865475))
f3 = Constant((0.0, 1.0))
def MaterialModel_3(X0, fiber, C0, C1, kappa_s):
    # 计算形变梯度
    FF = variable(grad(X0))
    J = det(FF)
    FF_bar = J**(-1/2)*FF
    CC_bar = FF_bar.T*FF_bar
    CC = FF.T*FF
    I1_bar = tr(CC_bar)
    m = fiber
    I4 = dot(m, CC*m)
    I4_bar = dot(m, CC_bar*m)
    # 本构参数
    Psi = 0.5*C0*(I1_bar-3) + C1*(exp(I4_bar-1) - I4_bar) + 0.5*kappa_s*(0.5*(J*J-1)-ln(J))
    dPsi = diff(Psi, FF)    
    dVs = TestFunction(Vs)
    L_hat = - inner(dPsi, grad(dVs))*dx
    return L_hat


if material == "n1":
    L_hat = MaterialModel_1(X0)
elif material == "n2":
    L_hat = MaterialModel_2(X0)
elif material == "f0":
    L_hat = MaterialModel_3(X0, f1, C0=C0, C1=0.0, kappa_s=kappa_s)
elif material == "f1":
    L_hat = MaterialModel_3(X0, f1, C0=C0, C1=C1, kappa_s=kappa_s)
elif material == "f2":
    L_hat = MaterialModel_3(X0, f2, C0=C0, C1=C1, kappa_s=kappa_s)
elif material == "f3":
    L_hat = MaterialModel_3(X0, f3, C0=C0, C1=C1, kappa_s=kappa_s)

L_hat += beta*inner(X_start-X0, TestFunction(Vs))*ds(1)

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




















