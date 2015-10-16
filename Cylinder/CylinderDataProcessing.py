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
# Main script to handle data processing
# **********************************************************

# Set up file IO
basedir = "Cylinder/";
Re = [52.3455, 61.5383, 75, 88.4617, 97.6545];
Xnew = np.array;
iter = 0;
for i in range(0,5):
    FileNameBase = basedir + "CylInitRe" + str(Re[i]);
    for j in range(0,6):
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
            X = np.zeros((data.nx*data.ny,30));
            X[0:data.nx*data.ny,iter] = Xnew;
        else:
            X[0:data.nx*data.ny,iter] = Xnew;
        iter = iter+1;
# POD
print "COMPUTING POD...\n"
Xmean = np.sum(X,axis=0)/np.size(X,1);
X = X - Xmean;
XTX = np.dot(np.transpose(X),X);
mu,v = np.linalg.eig(XTX);
M = 5;
PHI = np.zeros((np.size(X,0),M));
LAM = np.zeros(M);
for i in range(0,M):
    PHI[0:np.size(X,0),i] = 1/np.sqrt(mu[i])*np.dot(X,v[0:np.size(X,1),i]);
    LAM[i] = (1/M)*mu[i];
    contourPlot(data.X,data.Y,PHI[0:np.size(X,0),i]);
