from fenics import *
fluid_mesh = UnitSquareMesh(10,10)

#定义整个周期性边界条件
class PeriodicBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return bool(
            (x[0] < DOLFIN_EPS and on_boundary) or  # 左边界 x = 0
            (x[1] < DOLFIN_EPS and on_boundary)     # 下边界 y = 0
        )

    def map(self, x, y):
        # 检查右边界 x = 1 映射到 x = 0
        if near(x[0], 1.0):
            y[0] = x[0] - 1.0
            y[1] = x[1]
        # 检查上边界 y = 1 映射到 y = 0
        elif near(x[1], 1.0):
            y[0] = x[0]
            y[1] = x[1] - 1.0
        else:
            # 如果没有匹配到周期性条件，则保持原位置
            y[0] = x[0]
            y[1] = x[1]

V = FunctionSpace(fluid_mesh, "P", 1, constrained_domain= PeriodicBoundary())
V2 = FunctionSpace(fluid_mesh, "P", 1)

v = Function(V)

v.vector()[1] = 1.0
File("v.pvd") << v

v2 = Function(V2)

# v2.interpolate(v)
# v2.vector()[1] = 0.0
v2.vector()[1] = 1.0
File("v2.pvd") << v2

