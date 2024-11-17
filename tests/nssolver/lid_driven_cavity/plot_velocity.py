import pandas as pd
from ibfenics1.plot import plot_multiple_lines_1

def read_floats_from_file(file_path):
    float_list = []
    try:
        with open(file_path, 'r') as file:
            for line in file:
                # 去除行末的换行符并转换为浮点数
                float_list.append(float(line.strip()))
    except FileNotFoundError:
        print(f"文件未找到: {file_path}")
    except ValueError as e:
        print(f"文件中包含无法转换为浮点数的内容: {e}")
    return float_list

excel_file_s = [
    '/home/fenics/ibfenics/data/large/lid_driven_cavity_th_be/note/20241117-114939/results_x.xlsx',
    # '/home/fenics/ibfenics/data/large/lid_driven_cavity_chorin/note/20241117-120655/results_x.xlsx',
    # '/home/fenics/ibfenics/data/large/lid_driven_cavity_chorin/note/20241117-120906/results_x.xlsx',
    '/home/fenics/ibfenics/data/large/lid_driven_cavity_th_be/note/20241117-115842/results_x.xlsx',
    # '/home/fenics/ibfenics/data/large/lid_driven_cavity_th_be/note/20241117-120111/results_x.xlsx',
    ]

text_file_s = [
    # '/home/fenics/npuheart-1/build/u128.csv',
    '/home/fenics/npuheart-1/build/u256.csv',
    '/home/fenics/npuheart-1/build/u.csv',
    # '/home/fenics/npuheart-1/build/u_.csv',
]

x_list_s = []
p_list_s = []

for excel_file in excel_file_s:
    df = pd.read_excel(excel_file, sheet_name='y_list_x')
    x_list = df.iloc[:, 0].tolist()
    df = pd.read_excel(excel_file, sheet_name='u_list_x')
    p_list = df.iloc[:, 0].tolist()
    x_list_s.append(x_list)
    p_list_s.append(p_list)

for text_file in text_file_s:
    p_list = read_floats_from_file(text_file)
    p_list_s.append(p_list)
    x_list = [i/(len(p_list)-2)-1/(len(p_list)-2)/2 for i in range(len(p_list))]
    x_list_s.append(x_list)


plot_multiple_lines_1(
    x_list_s,
    p_list_s,
    legends=["fenics-th-100","fenics-th-50","fenics-th-10", "u128", "u256","u32","u_"],
)