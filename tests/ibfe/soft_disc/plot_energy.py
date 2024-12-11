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
    grep "total_energy(u0)" {log_path} > {output_path}/E.txt
    """
    os.system(command_extract)


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
            print(number)
    return numbers
            

base_path = "/home/fenics/ibfenics/data/large/"
log_path_s = [
    "oscilating_ring_base/0.01/20241210-175143",
    "oscilating_ring_base/0.005/20241210-174350",
    "oscilating_ring_be_sav/0.01/20241210-183913",
    "oscilating_ring_be_sav/0.01/20241210-192236",
    "oscilating_ring_be_sav/0.005/20241210-184707",
    "oscilating_ring_be_sav/0.005/20241210-190307"
]

dt_s = [
    0.01,
    0.005,
    0.01,
    0.01,
    0.005,
    0.005,
       ]

times = []
Es = []
T=1.2
for i in range(6):
    log_path = base_path + log_path_s[i]
    print(log_path)
    extract_data_from_log(log_path + "/info.log", log_path)
    E = extract_number(log_path+"/E.txt")
    print(E)
    time = [dt_s[i] * j for j in range(len(E))]
    times.append(time)
    Es.append(E)


plot_multiple_lines_1(
    times,
    Es,
    ylim=[1.0, 2],
    xlim=[0, T],
    # ylog=True,
    legends=[
        "IB-BE, $\Delta t=0.01$",
        "IB-BE, $\Delta t=0.005$",
        "SAV-IB-BE, $\Delta t=0.01,\\alpha=0$",
        "SAV-IB-BE, $\Delta t=0.01,\\alpha>0$",
        "SAV-IB-BE, $\Delta t=0.005,\\alpha=0$",
        "SAV-IB-BE, $\Delta t=0.005,\\alpha>0$",
        ],
    # markers=["o", "x", "s"],
    # labels=["bdf2-sav", "be"],
    title="",
    xlabel="",
    ylabel="",
    ncol=1,
)
