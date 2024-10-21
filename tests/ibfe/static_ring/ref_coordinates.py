import numpy as np
import warnings
from scipy.optimize import fsolve
from local_mesh import get_mesh
from fenics import *

class FiberForce(UserExpression):
    gamma = 0.0
    R = 0.25
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def eval(self, values, x):
        R, gamma = self.R, self.gamma
        s = self.fun(x[0], x[1])
        # print(s)
        # print(x[0], x[1])
        r = (-np.cos(s[0]/R), -np.sin(s[0]/R))
        values[0] = (1+s[1])/R*r[0]
        values[1] = (1+s[1])/R*r[1]

    def value_shape(self):
        return (2,)


    def fun(self, X0, X1):
        R, gamma = self.R, self.gamma
        def f(s2):
            a = X0 - 0.5
            b = X1 - 0.5
            c = R + s2
            d = R + s2 + gamma
            return 1 - (a / c)**2 - (b / d)**2

        initial_guess = -0.03
        
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")  # 记录所有警告
            # solution = fsolve(equation, 1.0)
            s2, info, ier, msg = fsolve(f, initial_guess, full_output=True, maxfev=100)
            # print(f"info: {info}, ier: {ier}, msg: {msg}")

            # 检查是否捕获到警告
            if w:
                for warning in w:
                    print(f"捕获到警告: {warning.message}")
        
        
        s1 = np.arccos((X0 - 0.5) / (R + s2)) * R
        s1 = 1.0 if np.isnan(s1) else s1  # arccos(?) = nan when ? > 1
        s1 = 2.0 * np.pi * R - s1 if X1 < 0.5 else s1
        return s1, s2



if __name__ == "__main__":
    solid_mesh = get_mesh(40)
    V = VectorFunctionSpace(solid_mesh, "CG", 1)
    velocity_inlet = FiberForce()
    u = Function(V)
    u.interpolate(velocity_inlet)
    File("velocity_inlet.pvd") << u




