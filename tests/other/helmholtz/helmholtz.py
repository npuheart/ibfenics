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



F = rho/dt*u*v*dx + mu*inner(grad(v), grad(u)) * dx
a = lhs(F)
L = rhs(F)
u = Function(V)

solve(a == L, u, bcs=bcu)
File("u.pvd") << u
u_list_x, x_list_x = extract_over_mid_x(u, [0.0,1.0],[0.0,1.0], n=100)

data_list = [{'u_list_x': u_list_x}, 
             {'x_list_x': x_list_x},]
write_excel_sheets(data_list, "a.xlsx", ["u_list_x",  "x_list_x" ])

