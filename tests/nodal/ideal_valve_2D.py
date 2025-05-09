import os
import swanlab
from loguru import logger
from fenics import *
from local_mesh import *
from ibfenics1.nssolver import IPCSSolver
from ibfenics1.io import (
    unique_filename,
    create_xdmf_file,
    TimeManager,
    write_paramters,
    write_excel,
    write_excel_sheets,
)




u0 = Function(V, name="velocity")
p0 = Function(Q, name="pressure")
f0 = Function(V, name="force")

time_manager = TimeManager(T, num_steps, 20)
bcu, bcp = calculate_fluid_boundary_conditions(V, Q)
fluid_solver = IPCSSolver(u0, p0, f0, dt, nv, bcp=bcp, bcu=bcu, rho=1.0, conv=True)
file_xlsx_name_x = unique_filename(
    os.path.basename(__file__), "note", "/results_x.xlsx"
)
file_xlsx_name_y = unique_filename(
    os.path.basename(__file__), "note", "/results_y.xlsx"
)
file_json_name = unique_filename(
    os.path.basename(__file__), "note", "/results.json"
)
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



u_list_x, v_list_x, p_list_x, y_list_x = extract_over_mid_x(u0, p0, [0.0,1.0],[0.0,1.0], n=100)
u_list_y, v_list_y, p_list_y, x_list_y = extract_over_mid_y(u0, p0, [0.0,1.0],[0.0,1.0], n=100)

data_list = [{'u_list_y': u_list_y}, 
             {'v_list_y': v_list_y},
             {'p_list_y': p_list_y}, 
             {'x_list_y': x_list_y},]
write_excel_sheets(data_list, file_xlsx_name_y, ["u_list_y", "v_list_y", "p_list_y", "x_list_y" ])

data_list = [{'u_list_x': u_list_x}, 
             {'v_list_x': v_list_x},
             {'p_list_x': p_list_x}, 
             {'y_list_x': y_list_x},]
write_excel_sheets(data_list, file_xlsx_name_x, ["u_list_x", "v_list_x", "p_list_x", "y_list_x" ])