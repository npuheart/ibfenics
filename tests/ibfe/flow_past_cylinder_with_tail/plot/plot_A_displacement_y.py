import pandas as pd
from ibfenics1.plot import plot_multiple_lines_1



def read_coordinate_x(filename):
    df = pd.read_csv(filename)
    x_values = df["avg(displacement (1))"]
    x_values = [data for data in x_values]
    return x_values

filename = "/home/fenics/ibfenics/tests/ibfe/flow_past_cylinder_with_tail/plot/b.csv"
a = read_coordinate_x(filename)

# def read_first_column(path_excel):
#     path_root = "/home/fenics/ibfenics/data/large/"
#     df = pd.read_excel(path_root + path_excel)
#     return df.iloc[:, 0].tolist()


# path_excel_s = [
#     "soft_disc_base/0.005/20241110-230541/volume.xlsx",
# ]
# volume_list_s = [read_first_column(path_excel) for path_excel in path_excel_s]
# volume_0 = volume_list_s[0][0]
# lost_percentage_s = [
#     [abs((volume - volume_0) / volume_0) * 100 for volume in volume_list]
#     for volume_list in volume_list_s
# ]
dt = 0.005*2
# dt = 0.005/6
time = [dt * i for i in range(len(a))]
times = [time]
plot_multiple_lines_1(
    times,
    [a],
    # ylim=[1e-10, 1e1],
    xlim=[0, 12],
    # ylog=True,
    # legends=["P3", "P4"],
    # labels=["$R_j^n$"],
    linestyles=["-"],
    legends=False,
    # markers=["o", "x", "s"],
    # labels=["bdf2-sav", "be"],
    title="",
    xlabel="time($s$)",
    ylabel="$y_A(m)$",
    ncol=1,
)
