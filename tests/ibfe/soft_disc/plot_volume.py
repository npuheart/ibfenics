import pandas as pd
from plot import plot_multiple_lines_1


def read_first_column(path_excel):
    path_root = "/home/fenics/ibfenics1/data/large/"
    df = pd.read_excel(path_root + path_excel)
    return df.iloc[:, 0].tolist()


path_excel_s = [
    "soft_disc_bdf2_sav/0.005/20241110-225536/volume.xlsx",
    "soft_disc_bdf2_sav/0.005/20241110-225536/volume.xlsx",
    "soft_disc_bdf2_sav/0.005/20241110-225536/volume.xlsx",
    "soft_disc_bdf2_sav/0.005/20241110-225536/volume.xlsx",
    "soft_disc_bdf2_sav/0.005/20241110-225536/volume.xlsx",
    "soft_disc_bdf2_sav/0.005/20241110-225536/volume.xlsx",
]
volume_list_s = [read_first_column(path_excel) for path_excel in path_excel_s]
volume_0 = volume_list_s[0][0]
lost_percentage_s = [
    [abs((volume - volume_0) / volume_0) * 100 for volume in volume_list]
    for volume_list in volume_list_s
]
dt = 0.005
time = [0.005 * i for i in range(len(volume_list_s[0]))]
times = [time for _ in range(len(lost_percentage_s))]
plot_multiple_lines_1(
    times,
    lost_percentage_s,
    ylim=[1e-10, 1e1],
    xlim=[0, 10],
    ylog=True,
    # legends=["P3", "P4"],
    labels=["bdf2-sav", "bdf2-dual", "be-sav", "be-dual", "bdf2", "be"],
    # labels=["bdf2-sav", "be"],
    title="Volume",
    xlabel="Time",
    ylabel="Volume",
    ncol=1,
)
