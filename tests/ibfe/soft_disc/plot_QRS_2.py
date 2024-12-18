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
    grep "150 - S " {log_path} > {output_path}/S.txt
    grep "151 - R " {log_path} > {output_path}/R.txt
    grep "152 - Q " {log_path} > {output_path}/Q.txt
    """
    os.system(command_extract)


figname = "1"
title="SAV-IB-BDF2 $\Delta t=0.00125$ $\\alpha=0.0$"
log_path = "/home/fenics/ibfenics/data/large/oscilating_ring_bdf2_sav/0.00125/20241216-175327"


# figname = "2"
# title="SAV-IB-BDF2 $\Delta t=0.00125$ $\\alpha>0.0$"
# log_path = "/home/fenics/ibfenics/data/large/oscilating_ring_bdf2_sav/0.00125/20241217-144409"

# figname = "3"
# title="SAV-IB-BDF2 $\Delta t=0.0025$ $\\alpha=0.0$"
# log_path = "/home/fenics/ibfenics/data/large/oscilating_ring_bdf2_sav/0.0025/20241216-175841"

# figname = "4"
# title="SAV-IB-BDF2 $\Delta t=0.0025$ $\\alpha>0.0$"
# log_path = "/home/fenics/ibfenics/data/large/oscilating_ring_bdf2_sav/0.0025/20241217-144414"

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
    ylim=[0.9, 1.2],
    xlim=[0, T],
    # ylog=True,
    legends=["$R^n$", "$q^n$", "$S^n$"],
    # markers=["o", "x", "s"],
    # labels=["bdf2-sav", "be"],
    title=title,loc="lower left",
    figname=figname,
    xlabel="",
    ylabel="",
    aspect=1,
    ncol=1,
)
