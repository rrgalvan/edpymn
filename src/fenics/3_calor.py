from fenics import *
import matplotlib.pylab as plt

# Build mesh from an external file
mesh_file = "m.xdmf"
mesh = Mesh()
with XDMFFile(mesh_file) as file:
    file.read(mesh)

do_plot = True
if do_plot:
    plot(mesh)
    plt.show()

Vh = FunctionSpace(mesh, 'P', 2) # Espacio de elementos finitos
u, ub = TrialFunction(Vh), TestFunction(Vh)  # Unknown and test function

# Datos
nu = 10  # Coeficiente de difusión
dt = 0.01  # Paso en tiempo
nIter = 100  # Nº iteraciones en tiempo
Tmax = dt*nIter
t = 0 # Instante inicial

# Fuerza externa, f
C0 = 1; C1 = 30
f = Expression("C0*exp(-C1*(pow(x[0],2) + pow(x[1],2))) * t/Tmax",
               C0=C0, C1=C1, Tmax=Tmax, t=t, degree=2);

u0 = interpolate(f, Vh)  # Condición inicial

vtkFile = File("calor.pvd")
u0.rename("u", "calor")  # Nombre que veremos en Paraview
vtkFile << (u0, t)

# Formulación variacional
a = u*ub*dx + dt*dot(grad(u),grad(ub))*dx
L = u0*ub*dx + f*ub*dx

# Bucle en tiempo (Euler Implícito)
u = Function(Vh)
u.rename("u", "calor")
for m in range(1, nIter+1):  # m = 1, 2,..., nIter
    # Actualizar el tiempo actual
    t = dt*m;
    f.t = t
    print(f"Iteraci'on de tiempo {m} (t={t})")

    solve(a==L, u)
    if do_plot:
        p = plot(u, cmap="coolwarm")
        plt.colorbar(p)
        plt.show()

    # Grabamos en fichero vtk
    vtkFile << (u, t)

    # Preparamos la siguiente etapa
    u0.assign(u)
