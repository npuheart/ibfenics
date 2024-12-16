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
    grep "total_energy(u0) " {log_path} > {output_path}/E.txt
    """
    # grep "total_energy(u0)" {log_path} > {output_path}/E.txt
    
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
    "soft_disc_bdf2/0.025/20241216-123254",
    "soft_disc_bdf2/0.0125/20241216-123600",
    "soft_disc_bdf2_sav/0.025/20241216-135653",
    "soft_disc_bdf2_sav/0.0125/20241216-124712",
]
dt_s = [
    0.0125,
    0.025,
    0.025,
    0.0125,
]
times = []
Es = []
T=12
for i in range(4):
    log_path = base_path + log_path_s[i]
    print(log_path)
    extract_data_from_log(log_path + "/info.log", log_path)
    E = extract_number(log_path+"/E.txt")
    print(E)
    time = [dt_s[i] * j for j in range(len(E))]
    times.append(time)
    Es.append(E)


# Es = [Es[0]]
# times = [times[0]]
plot_multiple_lines_1(
    times,
    Es,
    ylim=[0.0, 0.1],
    xlim=[0, T],
    # ylog=True,
    legends=[
        "IB-BDF2, $\Delta t=0.025$",
        "IB-BDF2, $\Delta t=0.0125$",
        "SAV-IB-BDF2, $\Delta t=0.025$",
        "SAV-IB-BDF2, $\Delta t=0.0125$",
        ],
    # markers=["o", "x", "s"],
    # labels=["bdf2-sav", "be"],
    linestyles = ['-.', '--','-',  ':'],
    title="",
    xlabel="",
    ylabel="",
    loc="upper right",
    ncol=1,
)
