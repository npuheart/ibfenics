

# double pressure_analytical(const double3& x) {
#     double p0      = 0.0;
#     double mu      = 1.0;
#     double R       = 0.25;
#     double w       = 0.0625;
#     double r       = std::sqrt((x.x - 0.5) * (x.x - 0.5) + (x.y - 0.5) * (x.y - 0.5)); // 计算到圆心的距离
#     double s2      = r - 0.25;
#     double result  = p0 - (mu / (w * R)) * s2 - 0.5 * (mu / (w * R)) * s2 * s2;
#     double result2 = p0 - (mu / R) - 0.5 * (mu * w / R);
#     result > 0 ? result = 0 : result = result;
#     result < result2 ? result = result2 : result = result;
#     return result;
# }

import numpy as np
from fenics import UserExpression

class PressureExact(UserExpression):
    #
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
    #
    def eval(self, values, x):
        values[0] = self.fun(x[0], x[1])
    #
    def value_shape(self):
        return ()
    #
    def fun(self, x, y):
        p0      = 0.0
        mu      = 1.0
        R       = 0.25
        w       = 0.0625
        r       = np.sqrt((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5))
        s2      = r - 0.25
        result  = p0 - (mu / (w * R)) * s2 - 0.5 * (mu / (w * R)) * s2 * s2
        result2 = p0 - (mu / R) - 0.5 * (mu * w / R)
        result = 0 if result > 0 else result
        result = result2 if result < result2 else result
        return result


# from fenics import UnitSquareMesh, FunctionSpace, Function, File
# from local_mesh import solid_mesh

# print(f"solid_mesh.id(): {solid_mesh.id()}")

# mesh = UnitSquareMesh(64,64)
# V = FunctionSpace(solid_mesh,"P",1)
# f = Function(V)


pressure_exact = PressureExact()

from fenics import assemble, dx, Expression, Function, interpolate, inner
import numpy as np


max_abs_function = lambda f: np.max(np.abs(f.vector()[:]))
min_abs_function = lambda f: np.min(np.abs(f.vector()[:]))

def calculate_error_p(p0):
    p1 = interpolate(pressure_exact, p0.function_space())
    p0_ = Function(p0.function_space())
    p1_ = Function(p1.function_space())
    p0_.vector()[:] = p0.vector()[:]
    p1_.vector()[:] = p1.vector()[:]
    p0_.vector()[:] = p0_.vector()[:] - assemble(p0*dx)
    p1_.vector()[:] = p1_.vector()[:] - assemble(p1*dx)
    norm_2 = np.sqrt(assemble((p0_ - p1_)**2*dx))
    norm_1 = assemble(abs(p0_ - p1_)*dx)
    p0_.vector()[:] = p0_.vector()[:] - p1_.vector()[:]
    norm_infinity = max_abs_function(p0_)
    return norm_1, norm_2, norm_infinity

def calculate_error_u(u):
    norm_2 = np.sqrt(assemble(inner(u,u)*2*dx))
    norm_1 = assemble((abs(u[0])+abs(u[1]))*dx)
    norm_infinity = max_abs_function(u)
    return norm_1, norm_2, norm_infinity


__all__ = ["calculate_error_u", "calculate_error_p"]
