import sys
import numpy as np
import pylab
import os.path
import shutil
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter, ScalarFormatter, MaxNLocator
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, axes, plot, xlabel, ylabel, title, grid, savefig, show
# User-defined auxiliary functions/structures
from contourPlot import contourPlot
from Bin2Dat import Bin2Dat
from IBPMData import IBPMData

# **********************************************************
# MatPlotLib formatting stuff
# **********************************************************

inches_per_pt = 1.0/72.27
width = 1150
height = 400
fig_size = [width*inches_per_pt,height*inches_per_pt]

params = {   
    #'axes.labelsize': 30,
    #'text.fontsize': 40,
    #'legend.fontsize': 20,
    'xtick.labelsize': 40,
    'ytick.labelsize': 40,             
    'figure.figsize':fig_size,
    #'figure.markersize': 50}
}
pylab.rcParams.update(params);
plt.ion();

# **********************************************************
# Subroutine for reshaping data snapshot
# **********************************************************

def reshapeSnapshot(Z,nx,ny):
    # Function to reshape a data snapshot to (ny x nx) array
    iter = 0;
    Znew = np.zeros((ny,nx));
    for i in range(0,nx):
        for j in range(0,ny):
            Znew[j,i] = Z[iter];
            iter = iter+1;
    return Znew;

# **********************************************************
# Main script to handle data processing
# **********************************************************

# Set up file IO
basedir = "Cylinder/";
Re = [52.3455, 61.5383, 75, 88.4617, 97.6545];
Xnew = np.array;
iter = 0;
runs = np.size(Re);
times = 7;
snapshots = runs*times;
for i in range(0,runs):
    FileNameBase = basedir + "CylInitRe" + str(Re[i]);
    for j in range(0,times):
        ind = 2*(j+1);
        if j<4:
            FileName = FileNameBase + "/ibpm" + str(0) + str(ind) + "000.bin";
        else:
            FileName = FileNameBase + "/ibpm" + str(ind) + "000.bin";
        # Scan in data
        print "PROCESSING " + FileName + "...\n"
        data = IBPMData();
        data = Bin2Dat(FileName);
        # Append data as column of POD matrix
        Xnew = np.reshape(data.Z,data.nx*data.ny,1);
        if ((i==0) & (j==0)):
            X = np.zeros((data.nx*data.ny,snapshots));
            X[0:data.nx*data.ny,iter] = Xnew;
        else:
            X[0:data.nx*data.ny,iter] = Xnew;
        iter = iter+1;
nx = data.nx; ny = data.ny;
sizeX = np.size(X,0);
# Extract POD modes
print "COMPUTING POD...\n";
Xmean = np.sum(X,axis=1)/snapshots;
X = X - np.transpose(np.tile(Xmean,[snapshots,1]));
XTX = np.dot(np.transpose(X),X);
mu,v = np.linalg.eig(XTX);
M = 5;
PHI = np.zeros((sizeX,M));
LAM = np.zeros(M);
for i in range(0,M):
    PHI[0:sizeX,i] = 1/np.sqrt(mu[i])*np.dot(X,v[0:snapshots,i]);
    LAM[i] = (1/M)*mu[i];
    #contourPlot(data.X,data.Y,reshapeSnapshot(PHI[0:sizeX,i],data.nx,data.ny));
# Project data onto POD modes
coeff = np.zeros((M,snapshots));
for i in range(0,snapshots):
    for j in range(0,M):
        coeff[j,i] = np.dot(np.transpose(PHI[0:sizeX,j]),X[0:sizeX,i]);
fig = plt.figure()
ax = plt.gca()
for i in range(0,runs):
    indStart = times*i;
    indEnd = times*(i+1);
    plt.figure();
    for j in range(0,M):
        plt.plot(np.linspace(1,times,times),coeff[j,indStart:indEnd],marker="o");

# Wait for user to manually end program
_ = raw_input("Press [enter] to continue.");
