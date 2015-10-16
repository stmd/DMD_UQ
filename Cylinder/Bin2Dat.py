import sys
import numpy as np
import pylab
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter, ScalarFormatter, MaxNLocator
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, axes, plot, xlabel, ylabel, title, grid, savefig, show

# **********************************************************
# Struct for handling data IO
# **********************************************************
class Data:
    def __init__(self):
        self.nx = 0;
        self.ny = 0;
        self.ngrid = 0;
        self.dx = 0;
        self.x0 = 0;
        self.y0 = 0;
        self.numPoints = 0;
        self.q = 0;
        self.omega = 0;
        self.FX = 0;
        self.FY = 0;
        self.timestep = 0;
        self.time = 0;

# **********************************************************
# Function to read IBPM binary data and convert to .dat file
# **********************************************************
def Bin2Dat(filename):
    # Open file
    fid = open(filename, "r")
    try:
        # Read grid info
        nx = np.fromfile(fid,dtype=np.uint32,count=1);
        ny = np.fromfile(fid,dtype=np.uint32,count=1);
        ngrid = np.fromfile(fid,dtype=np.uint32,count=1);
        dx = np.fromfile(fid,dtype=np.float64,count=1);
        x0 = np.fromfile(fid,dtype=np.float64,count=1);
        y0 = np.fromfile(fid,dtype=np.float64,count=1);
        numPoints = np.fromfile(fid,dtype=np.uint32,count=1);
        # Preallocations
        numFluxes = 2*nx*ny + nx + ny;
        q = np.zeros((ngrid,numFluxes));
        omega = np.zeros((ngrid,nx,ny));
        FX = np.zeros(numPoints);
        FY = np.zeros(numPoints);
        # Read flux
        print "Reading flux..."
        for lev in range(0,ngrid):
            for qind in range(0,numFluxes):
                q[lev,qind] = np.fromfile(fid,dtype=np.float64,count=1);
        # Read vorticity
        print "Reading vorticity..."
        for lev in range(0,ngrid):
            for i in range(1,nx):
                for j in range(1,ny):
                    omega[lev,i,j] = np.fromfile(fid,dtype=np.float64,count=1);
        # Read boundary force
        print "Reading boundary force vectors..."
        for i in range(0,numPoints):
            FX[i] = np.fromfile(fid,dtype=np.float64,count=1);
            FY[i] = np.fromfile(fid,dtype=np.float64,count=1);
        # Read timestep and time
        print "Reading timestep and time..."
        timestep = np.fromfile(fid,dtype=np.float64,count=1);
        time = np.fromfile(fid,dtype=np.float64,count=1);
        # Return struct containing information
        data = Data();
        data.nx = nx;
        data.ny = ny;
        data.ngrid = ngrid;
        data.dx = dx;
        data.x0 = x0;
        data.y0 = y0;
        data.numPoints = numPoints;
        data.q = q;
        data.omega = omega;
        data.FX = FX;
        data.FY = FY;
        data.timestep = timestep;
        data.time = time;
        
        return data;

    finally:
        fid.close()


# **********************************************************
# Main script to handle data processing
# **********************************************************

# Read file data
FileName = "Cylinder/CylInitRe52.3455/ibpm10000.bin";
data = Bin2Dat(FileName);
# Display contour plot
inches_per_pt = 1.0/72.27
width = 1150
height = 400
fig_size = [width*inches_per_pt,height*inches_per_pt]

params = {   #'axes.labelsize': 30,
             #'text.fontsize': 40,
             #'legend.fontsize': 20,
             'xtick.labelsize': 40,
             'ytick.labelsize': 40,             
             'figure.figsize':fig_size,
             #'figure.markersize': 50}
}
pylab.rcParams.update(params)

fig = plt.figure()
ax = plt.gca()
# Colorbar code ************
min, max = (0,100)
step = 1
# Using contourf to provide my colorbar info, then clearing the figure
ZZ = [[0,0],[0,0]]
levels = range(min,max+step,step)
CS3 = plt.contourf(ZZ, levels, cmap=cm.coolwarm)
plt.clf()
# Set up rectangular grid of (x,y) points
x0 = data.x0; y0 = data.y0;
nx = data.nx; ny = data.ny;
dx = data.dx;
omega = data.omega;
xf = x0 + (nx-1)*dx;
yf = y0 + (ny-1)*dx;
x = np.linspace(x0,xf,nx);
y = np.linspace(y0,yf,ny);
Z = np.squeeze(omega[0,0:nx,0:ny]);
XX,YY = np.meshgrid(x,y);
XX = np.transpose(XX);
YY = np.transpose(YY);
X = np.zeros(nx*ny); Y = np.zeros(nx*ny);
iter = 0;
for i in range(0,nx):
    for j in range(0,ny):
        X[iter] = XX[i,j];
        Y[iter] = YY[i,j];
        iter = iter+1;
plt.scatter(X,Y,c=Z,cmap=cm.coolwarm,s=5,lw=0,vmin=-5,vmax=5);
ax2 = plt.gca();
ax2.set_aspect('equal')
cbar = plt.colorbar(CS3, ticks=[0, 100])
cbar.ax.set_yticklabels(['-5', '5'])

ax.set_axisbelow(True)
plt.gcf().subplots_adjust(left=0.16)
plt.gcf().subplots_adjust(bottom=0.15)
ax.xaxis.labelpad = 20
plt.gcf().tight_layout()
#plt.grid(None)
pylab.savefig('2Dplot.eps',bbox_inches=0)
#plt.gca().tight_layout()
plt.draw()
plt.show()
