import pandas as pd
from ibfenics1.plot import plot_multiple_lines_1

# 读取名为 'Sheet1' 的工作表
excel_file_s = [
    '/home/fenics/ibfenics/data/large/lid_driven_cavity_ipcs/note/20241114-223051/results_x.xlsx',
    '/home/fenics/ibfenics/data/large/lid_driven_cavity_ipcs/note/20241114-223100/results_x.xlsx',
    '/home/fenics/ibfenics/data/large/lid_driven_cavity_ipcs/note/20241114-223108/results_x.xlsx',
    '/home/fenics/ibfenics/data/large/lid_driven_cavity_ipcs/note/20241114-223116/results_x.xlsx',
    '/home/fenics/ibfenics/data/large/lid_driven_cavity_ipcs/note/20241114-223128/results_x.xlsx',
    '/home/fenics/ibfenics/data/large/lid_driven_cavity_ipcs/note/20241114-223255/results_x.xlsx',
    # '/home/fenics/ibfenics/data/large/lid_driven_cavity_ipcs/note/20241114-222812/results_x.xlsx',
    ]

x_list_s = []
p_list_s = []
for excel_file in excel_file_s:
    df = pd.read_excel(excel_file, sheet_name='y_list_x')
    x_list = df.iloc[:, 0].tolist()
    df = pd.read_excel(excel_file, sheet_name='v_list_x')
    p_list = df.iloc[:, 0].tolist()
    x_list_s.append(x_list)
    p_list_s.append(p_list)

plot_multiple_lines_1(
    x_list_s,
    p_list_s,
    legends=["8","16","24","32","40","50"],
)