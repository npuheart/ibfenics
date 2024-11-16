import os
from loguru import logger
from fenics import *
from local_mesh import *
from ibfenics1.io import (
    unique_filename,
    create_xdmf_file,
    TimeManager,
    write_paramters,
    write_excel,
    write_excel_sheets,
)

v = TestFunction(V)
u = TrialFunction(V)
un = Function(V)
F = rho*(u-un)/dt*v*dx + mu*inner(grad(v), grad(u)) * dx
a = lhs(F)
L = rhs(F)



F2 = inner(grad(v), grad(u)) * dx + rho/dt*(grad(un)[0]+grad(un)[1])*v*dx
a2 = lhs(F2)
L2 = rhs(F2)





u = Function(V)
p = Function(V)
for i in range(1):
    solve(a == L, u, bcs=bcu)
    un.assign(u)
    solve(a2 == L2, p, bcs=bcu)


File("u.pvd") << p
u_list_x, x_list_x = extract_over_mid_x(p, [0.0,1.0],[0.0,1.0], n=100)

print(u_list_x)
data_list = [{'u_list_x': u_list_x}, 
             {'x_list_x': x_list_x},]
write_excel_sheets(data_list, "a.xlsx", ["u_list_x",  "x_list_x" ])

