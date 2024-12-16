import pandas as pd
from ibfenics1.plot import plot_multiple_lines_1

def extract_number(log_file):
    import re
    numbers = []
    with open(log_file, "r") as file:
        # 逐行读取文件内容，并将每行数据添加到列表中
        lines = file.readlines()
        for line in lines:
            matches = re.search(r'\d+\.\d+$', line)
            if matches:
                number = float(matches.group())
                numbers.append(number)
            else:
                number = None
            print(number)
    return numbers
            

def read_first_column(path_excel):
    path_root = "/home/fenics/ibfenics/data/large/"
    df = pd.read_excel(path_root + path_excel)
    return df.iloc[:, 0].tolist()


path_excel_s = [
    "soft_disc_bdf2/0.0125/20241215-170216/volume.xlsx",
    "soft_disc_base/0.0125/20241215-170753/volume.xlsx",
    # "soft_disc_bdf2_sav/0.0125/20241215-203002/volume.xlsx", # alpha = 0
    "soft_disc_bdf2_sav/0.0125/20241215-214005/volume.xlsx", # alpha = 10
    # "soft_disc_be_sav/0.0125/20241215-210735/volume.xlsx",# alpha = 0
    "soft_disc_be_sav/0.0125/20241215-213113/volume.xlsx",# alpha = 10
    "soft_disc_bdf2_sav/0.025/20241215-221047/volume.xlsx",# alpha = 10, dt = 0.025
    # "soft_disc_base/0.005/20241110-230541/volume.xlsx",
    # "soft_disc_bdf2_sav/0.005/20241110-225536/volume.xlsx",
    # "soft_disc_base/0.005/20241110-230541/volume.xlsx",
    # "soft_disc_bdf2_sav/0.005/20241110-225536/volume.xlsx",
    # "soft_disc_base/0.005/20241110-230541/volume.xlsx",
    # "soft_disc_bdf2_sav/0.005/20241110-225536/volume.xlsx",
]


import os
output_path = "/home/fenics/ibfenics/data/large/soft_disc_bdf2/0.025/20241215-222520"
command_extract = f"""
grep "volume" {output_path}/info.log > {output_path}/E.txt
"""
os.system(command_extract)



volume_list_s = [read_first_column(path_excel) for path_excel in path_excel_s]
volume_list_s.append(extract_number(f"{output_path}/E.txt"))
volume_0 = volume_list_s[0][0]
lost_percentage_s = [
    [abs((volume - volume_0) / volume_0) * 100 for volume in volume_list]
    for volume_list in volume_list_s
]


dt_s = [0.0125, 0.0125, 0.0125, 0.0125,  0.025, 0.025]
times = []
for j in range(len(dt_s)):
    time = [dt_s[j] * i for i in range(len(volume_list_s[j]))]
    times.append(time)

plot_multiple_lines_1(
    times,
    lost_percentage_s,
    ylim=[1e-10, 1e1],
    xlim=[0, 10],
    ylog=True,
    legends=[
        "IB-BDF2$\quad\quad\quad\Delta t=0.0125$", 
        "IB-BE$\;\quad\quad\quad\quad\Delta t=0.0125$", 
        "SAV-IB-BDF2$\;\;\Delta t=0.0125$", 
        "SAV-IB-BE$\quad\quad\Delta t=0.0125$", 
        "SAV-IB-BDF2$\;\;\Delta t=0.025$", 
        "IB-BDF2$\quad\quad\quad\Delta t=0.025$"],
    # labels=["bdf2-sav", "be"],
    title="",
    xlabel="Time ($s$)",
    ylabel="Area change (%)",
    ncol=1,
)
