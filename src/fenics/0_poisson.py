from fenics import *
import matplotlib.pyplot as plt

nIntervals = 76;  # Number of intervals on each side of the boundary
# mesh = UnitSquareMesh(nIntervals, nIntervals)
# mesh = RectangleMesh(Point(1, 1), Point(3, 7), nIntervals, nIntervals)
mesh = Mesh("/tmp/rectan.xml")
plot(mesh); plt.show()
Vh = FunctionSpace(mesh, 'P', 1)  # Define P1 finite element space
def boundary(x, on_boundary):  # Function returning true for all x on boundary
    return on_boundary
bc = DirichletBC(Vh, 0, boundary)  # Homog. DirichletBC on all the boundary
u, v = TrialFunction(Vh), TestFunction(Vh)  # Unknown and test function
f = Expression("-2*(10-4*x[0]+pow(x[0],2)-8*x[1]+pow(x[1],2))", degree=2) #Constant(4)  # Right hand side function
a = dot(grad(u), grad(v))*dx  # Bilinear form (LHS)
L = f*v*dx  # Linear form (RHS)
u = Function(Vh)  # Redefine u, where are going to store the solution in it
solve(a==L, u, bc)
plot(u, mode="warp"); plt.show()

print("máx u:", max(u.vector()))
A = assemble(a)
print(A.array())
