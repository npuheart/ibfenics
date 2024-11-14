import numpy as np
import warnings
from scipy.optimize import fsolve
from local_mesh import *
from fenics import *


class FiberForce(UserExpression):
    gamma = 0.0
    R = 0.25

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def eval(self, values, x):
        R, gamma = self.R, self.gamma
        s = self.fun(x[0], x[1])
        mu = 1.0
        omega = 0.0625
        # print(s)
        # print(x[0], x[1])
        r = (-np.cos(s[0] / R), -np.sin(s[0] / R))
        values[0] = mu / omega * (1 + s[1]) / R * r[0]
        values[1] = mu / omega * (1 + s[1]) / R * r[1]

    def value_shape(self):
        return (2,)

    def fun(self, X0, X1):
        R, gamma = self.R, self.gamma
        s2 = sqrt((X0-0.5)*(X0-0.5) + (X1-0.5)*(X1-0.5))-R
        s1 = np.arccos((X0 - 0.5) / (R + s2)) * R
        s1 = 1.0 if np.isnan(s1) else s1  # arccos(?) = nan when ? > 1
        s1 = 2.0 * np.pi * R - s1 if X1 < 0.5 else s1
        return s1, s2


if __name__ == "__main__":
    V = VectorFunctionSpace(solid_mesh, "CG", 1)
    velocity_inlet = FiberForce()
    u = Function(V)
    u.interpolate(velocity_inlet)
    File("velocity_inlet.pvd") << u
