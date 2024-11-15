import pandas as pd


def read_coordinate_x(filename):
    df = pd.read_csv(filename)
    x_values = df["Points:0"]
    x_values = [data for data in x_values]
    return x_values


def read_pressure(filename):
    df = pd.read_csv(filename)
    pressure_values = df["pressure"]
    pressure_values = [data for data in pressure_values]
    return pressure_values


filenames = [
    "/home/fenics/ibfenics/tests/ibfe/static_ring/data_tmp/0.2.csv",
    "/home/fenics/ibfenics/tests/ibfe/static_ring/data_tmp/0.3.csv",
    "/home/fenics/ibfenics/tests/ibfe/static_ring/data_tmp/0.4.csv",
    "/home/fenics/ibfenics/tests/ibfe/static_ring/data_tmp/0.5.csv",
]


coordinate_x = read_coordinate_x(filenames[0])

pressures = []
for filename in filenames:
    pressure = read_pressure(filename)
    pressures.append(pressure)


coordinate_x_s = [coordinate_x for i in range(len(pressures))]


if __name__ == "__main__":
    from ibfenics1.plot import plot_multiple_lines_1

    plot_multiple_lines_1(
        coordinate_x_s,
        pressures,
        xlim=[0, 1],
        ylim=[-1, 3.5],
        legends=["0.2", "0.3", "0.4", "0.5"],
    )
    
__all__ = ["coordinate_x_s", "pressures"]
