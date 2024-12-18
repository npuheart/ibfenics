import math
def calculate_convergence_order(Eh1, Eh2, h1, h2):
    """
    计算收敛阶 p。
    
    参数:
    Eh1 (float): 对应网格 h1 的误差。
    Eh2 (float): 对应网格 h2 的误差。
    h1 (float): 第一个网格的步长。
    h2 (float): 第二个网格的步长。
    
    返回:
    float: 计算出的收敛阶 p。
    """
    # 确保输入的h1, h2不能为零
    if h1 <= 0 or h2 <= 0:
        raise ValueError("网格步长 h1 和 h2 必须大于零")
    
    # 计算收敛阶 p
    p = math.log(Eh1 / Eh2) / math.log(h1 / h2)
    return p


# BDF2
print(calculate_convergence_order(1.68e-02, 6.13e-02, 1/64, 1/32))
print(calculate_convergence_order(7.68e-03, 1.68e-02, 1/96, 1/64))
print(calculate_convergence_order(4.37e-03, 7.68e-03, 1/128, 1/96))

print(calculate_convergence_order( 4.39e-07,8.30e-07, 1/64, 1/32))
print(calculate_convergence_order(2.94e-07, 4.39e-07, 1/96, 1/64))
print(calculate_convergence_order( 2.22e-07,2.94e-07,  1/128, 1/96))

# BE

print(calculate_convergence_order(4.65e-07,8.79e-07, 1/64, 1/32))
print(calculate_convergence_order(3.11e-07,4.65e-07, 1/96, 1/64))
print(calculate_convergence_order(2.34e-07,3.11e-07,  1/128, 1/96))




