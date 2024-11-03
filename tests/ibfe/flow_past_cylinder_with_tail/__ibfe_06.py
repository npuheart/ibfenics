from ibfenics import TaylorHoodSolver, Interaction
from dolfin import *
from mshr import *
from datetime import datetime
import os
import numpy as np

TaylorHoodSolver = TaylorHoodSolver.TaylorHoodSolver


class IU:
    g = 0.001  # 三个基本单位：质量、长度、时间
    s = 1
    cm = 0.01
    m = 100 * cm
    # 基本单位(米、千克)
    kg = 1000 * g
    #
    dyn = g * cm / s / s
    # 力(达因、牛顿)
    N = kg * m / s / s
    #
    Pa = N / m / m
    # / 压强(帕斯卡、毫米汞柱)
    kPa = 1000 * Pa
    mmHg = 133.3223684 * Pa
    #
    J = N * m
    # 能量(焦耳、卡)
    cal = 4.1868 * J
    #
    W = J / s
    # 功率(瓦)


class VelocityInlet(UserExpression):
    def __init__(self, H, **kwargs):
        self.t = 0
        self.H = H
        super().__init__(**kwargs)

    def eval(self, values, x):
        if self.t >= 2.0 * IU.s:
            scale = 1.0
        else:
            scale = 0.5 * (1.0 - cos(0.5 * np.pi * self.t))

        values[0] = 1.5 * U / 0.5 / H / 0.5 / H * x[1] * (H - x[1]) * scale
        values[1] = 0.0

    def value_shape(self):
        return (2,)


# output boundaries
class Boundary_1(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] < DOLFIN_EPS_LARGE and on_boundary


class Boundary_4(SubDomain):
    def inside(self, x, on_boundary):
        return x[0] > L * (1.0 - DOLFIN_EPS_LARGE) and on_boundary


class Boundary_2(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] < DOLFIN_EPS_LARGE and on_boundary


class Boundary_3(SubDomain):
    def inside(self, x, on_boundary):
        return x[1] > H * (1.0 - DOLFIN_EPS_LARGE) and on_boundary


marker_bg_left = 1
marker_bg_right = 4
marker_bg_upper = 3
marker_bg_lower = 2

H = 0.41 * IU.m
L = 2.2 * IU.m
rho_f = 1000 * IU.kg / IU.m / IU.m / IU.m
rho_s = 1000 * IU.kg / IU.m / IU.m / IU.m
# nv_f = 0.001   *IU.Pa*IU.s
nv_f = 0.001 * IU.Pa * IU.s
nv_s = 0.4
E_s = 5.6e6 * IU.Pa
U = 2 * IU.m / IU.s
Ma = 0.22
dh = 0.001 * IU.m
dt = 1.0e-6 * IU.s

# Define fluid solver
beta = 5e11
rho = 1.0
T = 10.0

# stabliization parameter
alpha = 1.0 * dt
stab = False

current_file_name = os.path.basename(__file__)
note = os.path.splitext(current_file_name)[0]
file_id = "data/" + note + "-" + datetime.now().strftime("%Y%m%d-%H%M%S")

order_velocity = 2
order_pressure = 1
order_displacement = 1

solid_mesh = Mesh()
with XDMFFile("mesh/bar_behind_cylinder_40.xdmf") as xdmf:
    xdmf.read(solid_mesh)

bdry = MeshFunction("size_t", solid_mesh, "mesh/bar_behind_cylinder_40_boundaries.xml")
domains = MeshFunction("size_t", solid_mesh, "mesh/bar_behind_cylinder_40_domains.xml")
dx = Measure("dx")(subdomain_data=domains)

orders = [order_velocity, order_pressure, order_displacement]
points = [Point(0, 0), Point(L, H)]
seperations = [256, int(256 * H / L)]
interaction = Interaction(points, seperations, solid_mesh, orders)

fluid_mesh = interaction.fluid_mesh
ib_mesh = interaction.ib_mesh
ib_interpolation = interaction.ib_interpolation
Vs = interaction.Vs
Vf = interaction.Vf
Vp = interaction.Vp

boundary_bg = MeshFunction("size_t", fluid_mesh, 1)
Boundary_1().mark(boundary_bg, marker_bg_left)
Boundary_2().mark(boundary_bg, marker_bg_lower)
Boundary_3().mark(boundary_bg, marker_bg_upper)
Boundary_4().mark(boundary_bg, marker_bg_right)

print("solid_mesh.hmax() ", solid_mesh.hmax())
print("solid_mesh.hmin() ", solid_mesh.hmin())
print("fluid_mesh.hmax() ", fluid_mesh.hmax())
print("fluid_mesh.hmin() ", fluid_mesh.hmin())
print("ratio(>2) = ", fluid_mesh.hmin() / order_velocity / solid_mesh.hmax())

# Define trial and test functions for solid
us = TrialFunction(Vs)
vs = TestFunction(Vs)

# Create functions for fluid
u0 = Function(Vf, name="velocity")
p0 = Function(Vp, name="pressure")
f = Function(Vf, name="force")

# Create functions for solid
velocity = Function(Vs, name="velocity")
disp = Function(Vs, name="displacement")
force = Function(Vs, name="force")

disp.interpolate(Expression(("x[0]", "x[1]"), degree=2))
ib_interpolation.evaluate_current_points(disp._cpp_object)

# Define interpolation object and fluid solver object
navier_stokes_solver = TaylorHoodSolver(
    u0,
    p0,
    f,
    dt,
    nv_f,
    stab=stab,
    alpha=alpha,
    bdry=boundary_bg,
    bc_Neumann=[marker_bg_right],
)
W = navier_stokes_solver.W

# # Define boundary conditions
# noslip = DirichletBC(W.sub(0), (0, 0), "near(x[0], 2.2) || near(x[0], 0.0) || near(x[1], 0.0)")
# upflow = DirichletBC(W.sub(0), (1, 0), "near(x[1], 0.41)")
# pinpoint = DirichletBC(W.sub(1), 0, "near(x[0], 0.0) && near(x[1], 0.0)", "pointwise")
# bcu = [noslip, upflow]
# bcp = [pinpoint]
# Define boundary conditions
velocity_inlet = VelocityInlet(H)
bcu_1 = DirichletBC(W.sub(0), velocity_inlet, boundary_bg, marker_bg_left)
bcu_2 = DirichletBC(W.sub(0), Constant((0, 0)), boundary_bg, marker_bg_upper)
bcu_3 = DirichletBC(W.sub(0), Constant((0, 0)), boundary_bg, marker_bg_lower)
bcp_2 = DirichletBC(W.sub(1), Constant(0), boundary_bg, marker_bg_right)
bcus = [bcu_1, bcu_2, bcu_3]
bcps = [bcp_2]


# Define variational problem for solid
F = grad(disp)
EE = 0.5 * (F.T * F - Identity(2))
lambda_lame = nv_s * E_s / ((1.0 + nv_s) * (1.0 - 2.0 * nv_s))
mu_lame = E_s / (2.0 * (1.0 + nv_s))
SS = lambda_lame * tr(EE) * Identity(2) + 2.0 * mu_lame * EE
P = SS * F.T
F2 = inner(P, grad(vs)) * dx + inner(us, vs) * dx

# 固定
x_start = SpatialCoordinate(solid_mesh)
F2 += beta * inner(disp - x_start, vs) * dx(1)

a2 = lhs(F2)
L2 = rhs(F2)
A2 = assemble(a2)

# Output Directory name
directory = file_id + "-" + str(nv_s) + "-" + str(dt)
if stab:
    directory = directory + "-" + str(alpha)

print("output directory : ", directory)

# Create files for storing solution
file_solid = XDMFFile(directory + "/solid.xdmf")
file_solid.parameters["rewrite_function_mesh"] = False
file_solid.parameters["functions_share_mesh"] = True
file_solid.parameters["flush_output"] = True

file_fluid = XDMFFile(directory + "/fluid.xdmf")
file_fluid.parameters["rewrite_function_mesh"] = False
file_fluid.parameters["functions_share_mesh"] = True
file_fluid.parameters["flush_output"] = True


def elastic_energy(disp):
    F = grad(disp)
    C = F * F.T
    Ic = tr(C)
    J = det(F)
    return assemble(0.5 * nv_s * (Ic - 2 - 2 * ln(J)) * dx)


t = dt
En = 0.0  # elastic energy
num_steps = int(T / dt)
for n in range(1, num_steps + 1):
    velocity_inlet.t = t
    En = elastic_energy(disp)
    print(En)
    if np.isnan(En):
        print("The simulation is blowing up!")
        break
    print("energy: ", En)
    # step 1. calculate velocity and pressure
    # 计算流体的速度和压力
    navier_stokes_solver.update(u0, p0)
    u1, p1 = navier_stokes_solver.solve(bcus, bcps)
    u0.assign(u1)
    p0.assign(p1)
    # step 2. interpolate velocity from fluid to solid
    ib_interpolation.fluid_to_solid(u0._cpp_object, velocity._cpp_object)
    # step 3. calculate disp for solid and update current gauss points and dof points
    disp.vector()[:] = velocity.vector()[:] * dt + disp.vector()[:]
    ib_interpolation.evaluate_current_points(disp._cpp_object)
    # step 4. calculate body force.
    b2 = assemble(L2)
    solve(A2, force.vector(), b2)
    # step 5. interpolate force from solid to fluid
    ib_interpolation.solid_to_fluid(f._cpp_object, force._cpp_object)
    # step 6. update variables and save to file.
    file_fluid.write(u0, t)
    file_fluid.write(p0, t)
    file_fluid.write(f, t)
    file_solid.write(disp, t)
    file_solid.write(force, t)
    file_solid.write(velocity, t)
    t = n * dt
    print(t)
