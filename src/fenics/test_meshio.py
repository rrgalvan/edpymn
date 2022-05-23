from fenics import *
import matplotlib.pylab as plt

def convert_mesh_to_xdmf(input_file, output_file):
    import numpy as np
    import meshio
    msh = meshio.read(input_mesh_file)
    clls = np.vstack((
        cell.data for cell in msh.cells if cell.type == "triangle"
    ))  # only write 2D cells
    meshio.xdmf.write(output_mesh_file,
                      meshio.Mesh(msh.points, cells = {"triangle": clls}))

input_mesh_file = "/tmp/m.msh"
output_mesh_file = "/tmp/m.xdmf"
convert_mesh_to_xdmf(input_mesh_file, output_mesh_file)
mesh = Mesh()
with XDMFFile(output_mesh_file) as file:
    file.read(mesh)
plot(mesh)
plt.show()
