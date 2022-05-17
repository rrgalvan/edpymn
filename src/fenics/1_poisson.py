from fenics import *
import matplotlib.pyplot as plt

nIntervals = 40;  # Number of intervals on each side of the boundary
mesh = UnitSquareMesh(nIntervals, nIntervals)
do_plot = False
if do_plot:
    plot(mesh); plt.show()

Vh = FunctionSpace(mesh, 'P', 1)  # Define P1 finite element space
def boundary(x, on_boundary):  # Function returning true for all x on boundary
    return on_boundary
bc = DirichletBC(Vh, 0, boundary)  # Homog. DirichletBC on all the boundary
u, v = TrialFunction(Vh), TestFunction(Vh)  # Unknown and test function
f = Expression("2*x[0]*(1-x[0]) + 2*x[1]*(1-x[1])", degree=2)  # Right hand side function
a = dot(grad(u), grad(v))*dx  # Bilinear form (LHS)
L = f*v*dx  # Linear form (RHS)
u = Function(Vh)  # Redefine u, where are going to store the solution in it

solve(a==L, u, bc)

if do_plot:
    plot(u, mode="warp"); plt.show()

exactSol = Expression("x[0]*(1-x[0])*x[1]*(1-x[1])", degree=2)  # Exact soluton for previous data
errL2 = errornorm(exactSol, u, 'L2')
errH1 = errornorm(exactSol, u, 'H1')
print("Error L2:", errL2)
print("Error H1:", errH1)
