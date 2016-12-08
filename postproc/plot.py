import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
from vtk import *
from vtk.util.numpy_support import vtk_to_numpy

class MidpointNormalize(colors.Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

# load input data
reader = vtk.vtkXMLStructuredGridReader()
reader.SetFileName("output/solution82.vts")
reader.Update()

# Grab the scalar data
u_vtk = reader.GetOutput().GetPointData().GetArray("x mean velocity")
v_vtk = reader.GetOutput().GetPointData().GetArray("y mean velocity")

# Transform into numpy arrays
u = vtk_to_numpy(u_vtk)
v = vtk_to_numpy(v_vtk)
print("Size: "+str(u.shape[0]));
my = 100;
mx = 100;
u = np.reshape(u,[my,mx])
v = np.reshape(v,[my,mx])

speed = np.sqrt(u*u + v*v)

x = np.linspace(0,30.0,mx)
y = np.linspace(0,2.0,my)
X, Y = np.meshgrid(x,y)

# Plot the data
fig, ax = plt.subplots()
cf = ax.pcolormesh(X,Y,u, norm=MidpointNormalize(midpoint=0.),
                   cmap='RdBu_r')
strm = ax.streamplot(X,Y,u,v, color='k') #color = speed, cmap=plt.cm.autumn)
ax.set_xlim([np.min(x),np.max(x)])
ax.set_ylim([np.min(y),np.max(y)])
fig.set_size_inches(18,6,forward=True)

plt.show()

