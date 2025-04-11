from ibfenics1.cpp import IBMesh3D, IBInterpolation3D
from fenics import *
from mshr import Circle, generate_mesh





n_mesh = 32
points = [Point(0, 0, 0), Point(1, 2, 3)]
seperations = [n_mesh, n_mesh, n_mesh]
order_velocity = 2



ib_mesh = IBMesh3D(points, seperations, order_velocity)
fluid_mesh = ib_mesh.mesh()
solid_mesh = generate_mesh(Circle(Point(0.2, 0.2), 0.05), 64)
Vf = VectorFunctionSpace(fluid_mesh, "P", order_velocity)
uf = interpolate(Expression(("x[0]", "x[1]", "x[2]"), degree=2), Vf)


IBInterpolation3D(ib_mesh, solid_mesh)




ib_mesh.build_map(uf._cpp_object)



uf = interpolate(Expression(("x[0]*x[0]", "x[1]*x[1]", "x[2]*x[2]"), degree=2), Vf)
a = ib_mesh.evaluate(Point(0.51,0.52,0.53), uf._cpp_object)
print(a[0]-0.51*0.51)









