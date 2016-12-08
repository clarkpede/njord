import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib as mpl
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
reader.SetFileName("output/solution30.vts")
reader.Update()

# Grab the scalar data
data = reader.GetOutput()
dim = data.GetDimensions()
mx = dim[0]
my = dim[1]

# Transform into numpy arrays
u = vtk_to_numpy(data.GetPointData().GetArray('x mean velocity'))
v = vtk_to_numpy(data.GetPointData().GetArray('y mean velocity'))
u = np.reshape(u,[mx,my],order='F')
v = np.reshape(v,[mx,my],order='F')

speed = np.sqrt(u*u + v*v)

x = np.linspace(0,30.0,mx)
y = np.linspace(0,2.0,my)
Y, X = np.meshgrid(y,x)

# Plot the data
fig, ax = plt.subplots(2,1)
pcm1 = ax[0].pcolormesh(X,Y,u, norm=MidpointNormalize(midpoint=0.),
                       cmap='RdBu_r')
#strm = ax.streamplot(X,Y,u,v, color='k') #color = speed, cmap=plt.cm.autumn)
cbar1 = fig.colorbar(pcm1, ax = ax[0], orientation='vertical')
ax[0].set_xlim([np.min(x),np.max(x)])
ax[0].set_ylim([np.min(y),np.max(y)])


pcm2 = ax[1].pcolormesh(X,Y,v, norm=MidpointNormalize(midpoint=0.),
                       cmap='RdBu_r')
#strm = ax.streamplot(X,Y,u,v, color='k') #color = speed, cmap=plt.cm.autumn)
cbar2 = fig.colorbar(pcm2,ax = ax[1], orientation='vertical')
ax[1].set_xlim([np.min(x),np.max(x)])
ax[1].set_ylim([np.min(y),np.max(y)])

fig.set_size_inches(18,6,forward=True)
plt.show()

