from fenics import *
import matplotlib.pyplot as plt


nIntervals = 40;  # Number of intervals on each side of the boundary
mesh = UnitSquareMesh(nIntervals, nIntervals)
Vh = FunctionSpace(mesh, 'P', 1)

def boundary(x, on_boundary):  # Function returning true for all x on boundary
    return on_boundary

bc = DirichletBC(Vh, 0, boundary)  # Homog. DirichletBC on all the boundary
u, v = TrialFunction(Vh), TestFunction(Vh)  # Unknown and test function
f = Constant(4)  # Right hand side function
a = dot(grad(u), grad(v))*dx  # Bilinear form (LHS)
L = f*v*dx  # Linear form (RHS)

# Montar el sistema de ecuaciones
A = assemble(a)
b = assemble(L)
bc.apply(A,b)

# Otra posibiliad más sencilla
A, b = assemble_system(a, L, bc)

# Resolver
u = Function(Vh)
solve(A, u.vector(), b, "gmres", "ilu")

# Mostrar la gráfica
p = plot(u, cmap="coolwarm")  # Plot using a colormap https://matplotlib.org/3.5.0/tutorials/colors/colormaps.html
plt.colorbar(p)
plt.show()
