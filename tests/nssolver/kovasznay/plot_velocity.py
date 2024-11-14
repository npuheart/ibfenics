import pandas as pd
from ibfenics1.plot import plot_multiple_lines_1

# 读取名为 'Sheet1' 的工作表
excel_file_s = [
    '/home/fenics/ibfenics/data/large/kovasznay_th/note/20241114-172928/results.xlsx',
    '/home/fenics/ibfenics/data/large/kovasznay_ipcs/note/20241114-172919/results.xlsx']

x_list_s = []
p_list_s = []
for excel_file in excel_file_s:
    df = pd.read_excel(excel_file, sheet_name='x')
    x_list = df.iloc[:, 0].tolist()
    df = pd.read_excel(excel_file, sheet_name='p')
    p_list = df.iloc[:, 0].tolist()
    x_list_s.append(x_list)
    p_list_s.append(p_list)

plot_multiple_lines_1(
    x_list_s,
    p_list_s,
    legends=["TH", "IPCS"],
)