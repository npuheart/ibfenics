
from post_processing import *




pressure_exact = PressureExact()



from fenics import *


mesh = UnitSquareMesh(100, 100)
V = FunctionSpace(mesh, "CG", 1)
p = interpolate(pressure_exact, V)



def evaluate_along_line_y(p0, xlim=[0,1], line_y=0.5, n=10):
    p_list = []
    x_list = []
    for i in range(n):
        xx = xlim[0] + (xlim[1]-xlim[0])*i/(n-1)
        yy = line_y
        x_list.append(xx)
        p_list.append(p0(xx,yy)+4.125)
    
    return p_list, x_list 


line_y_s = [0.2,0.3,0.4,0.5]
p_list_s = []
x_list_s = []
for line_y in line_y_s:
    p_list, x_list = evaluate_along_line_y(p, [0,1], line_y, 100)
    p_list_s.append(p_list)
    x_list_s.append(x_list)    

from ibfenics1.plot import plot_multiple_lines_1

from plot_pressure import *

plot_multiple_lines_1(
    x_list_s + coordinate_x_s,
    p_list_s + pressures,
    xlim=[0, 1.1],
    # ylim=[-1, 3.5],
    linestyles=["solid","solid","solid","solid","dashed","dashed","dashed","dashed"],
    legends   =["0.2(exact)", "0.3(exact)", "0.4(exact)", "0.5(exact)","0.2", "0.3", "0.4", "0.5"],
    colors = ["#FABB6E",  "#ADDB88","#FAC7B3","#CEDFEF","#FC8002","#369F2D","#EE4431","#1663A9"])

