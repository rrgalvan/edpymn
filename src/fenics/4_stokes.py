from dolfin import *
import matplotlib.pyplot as plt

# Load mesh
# mesh = UnitSquareMesh.create(16, 16, CellType.Type.triangle)
mesh = UnitSquareMesh(16,16)
plot(mesh)
plt.show()

# Build function space
P2 = VectorElement("Lagrange", mesh.ufl_cell(), 2)
P1 = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
W = FunctionSpace(mesh, P2*P1)

# Next, we define the boundary conditions. ::

# Boundaries
def top(x, on_boundary):
    return x[1] > 1.0 - DOLFIN_EPS
def non_top(x, on_boundary):
    return on_boundary and x[1] <= 1.0 - DOLFIN_EPS

# No-slip boundary condition for velocity
noslip = Constant((0.0, 0.0))
bc0 = DirichletBC(W.sub(0), noslip, non_top)

# Inflow boundary condition for velocity
traction = Expression(("sin(x[1]*pi)", "0.0"), degree=2)
bc1 = DirichletBC(W.sub(0), traction, top)

# Collect boundary conditions
bcs = [bc0, bc1]

# The bilinear and linear forms corresponding to the weak mixed
# formulation of the Stokes equations are defined as follows: ::

# Define variational problem
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)
f = Constant((0.0, 0.0))
a = inner(grad(u), grad(v))*dx - div(v)*p*dx + q*div(u)*dx
L = inner(f, v)*dx

# Form for use in constructing preconditioner matrix
b = inner(grad(u), grad(v))*dx + p*q*dx

# Assemble system
A, b = assemble_system(a, L, bcs)

# Solve
U = Function(W)
solve(A, U.vector(), b)

# Finally, we can play with the result in different ways: ::

# Get sub-functions
u, p = U.split()

plot(u, title="Velocity")
plt.show()

plot(p, title="Pressure")
plt.show()


# Save solution in VTK format
ufile_pvd = File("velocity.pvd")
ufile_pvd << u
pfile_pvd = File("pressure.pvd")
pfile_pvd << p

import sys; sys.exit()
