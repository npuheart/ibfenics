import os, re
import pandas as pd
from ibfenics1.plot import plot_multiple_lines_1


def read_first_column(path_excel):
    path_root = "/home/fenics/ibfenics/data/large/"
    df = pd.read_excel(path_root + path_excel)
    return df.iloc[:, 0].tolist()


# 从日志提取数据
def extract_data_from_log(log_path, output_path):
    command_extract = f"""
    grep "153 - S " {log_path} > {output_path}/S.txt
    grep "154 - R " {log_path} > {output_path}/R.txt
    grep "155 - Q " {log_path} > {output_path}/Q.txt
    """
    os.system(command_extract)

# log_path = "/home/fenics/ibfenics/data/large/oscilating_ring_be_sav/0.01/20241210-041447/"

# title="SAV-IB-BE $\Delta t=0.005$ $\\alpha>0.0$"
# log_path = "/home/fenics/ibfenics/data/large/oscilating_ring_be_sav/0.005/20241210-190307"

# title="SAV-IB-BE $\Delta t=0.005$ $\\alpha=0.0$"
# log_path = "/home/fenics/ibfenics/data/large/oscilating_ring_be_sav/0.005/20241210-184707"

# title="SAV-IB-BE $\Delta t=0.01$ $\\alpha>0.0$"
# log_path = "/home/fenics/ibfenics/data/large/oscilating_ring_be_sav/0.01/20241210-192236"

title="SAV-IB-BE $\Delta t=0.01$ $\\alpha=0.0$"
log_path = "/home/fenics/ibfenics/data/large/oscilating_ring_be_sav/0.01/20241210-183913"

extract_data_from_log(log_path + "/info.log", log_path)


def extract_number(log_file):
    numbers = []
    with open(log_file, "r") as file:
        # 逐行读取文件内容，并将每行数据添加到列表中
        lines = file.readlines()
        for line in lines:
            matches = re.search(r'\d+\.\d+$', line)
            if matches:
                number = float(matches.group())
            else:
                number = None
            numbers.append(number)
            # print(number)
    return numbers
            

S = extract_number(log_path+"/S.txt")
Q = extract_number(log_path+"/Q.txt")
R = extract_number(log_path+"/R.txt")
print(S)
print(Q)
print(R)

T = 1.2
dt = T/len(S)
time = [dt * i for i in range(len(S))]
plot_multiple_lines_1(
    [time,time,time],
    [R, Q, S],
    ylim=[0.88, 1.18],
    xlim=[0, T],
    # ylog=True,
    legends=["$R^n$", "$q^n$", "$S^n$"],
    # markers=["o", "x", "s"],
    # labels=["bdf2-sav", "be"],
    title="SAV-IB-BE $\Delta t=0.01$ $\\alpha=0$",
    # title="SAV-IB-BE $\Delta t=0.01$ $\\alpha=10$",
    xlabel="",
    ylabel="",
    ncol=1,
)
