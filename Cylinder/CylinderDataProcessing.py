from __future__ import division
import sys
import numpy as np
import pylab
import os.path
import shutil
import time
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter, ScalarFormatter, MaxNLocator
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, axes, plot, xlabel, ylabel, title, grid, savefig, show
# User-defined auxiliary functions/structures
from contourPlot import contourPlot
from Bin2Dat import Bin2Dat
from readPltData import readPltData
from IBPMData import IBPMData
from LegendreProjection import Legendre, LegendreProjection

# **********************************************************
# MatPlotLib formatting stuff
# **********************************************************

inches_per_pt = 1.0/72.27
width = 600
height = 200
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
# Read/load data
# **********************************************************

# Parameters for runs
basedir = "Cylinder/";
Re = [52.3455, 61.5383, 75, 88.4617, 97.6545];
endPts = [50,100];
IC = [58000, 28000, 30000, 26000, 28000];
Xnew = np.array;
iter = 0;
runs = np.size(Re);
times = 100;
samp = 10;
DT = 0.02*samp;
snapshots = runs*times;
LoadRead = "READ";
# Either read data, or load from file
if LoadRead == "LOAD":
    for i in range(0,runs):
        # Get filename
        FileNameBase = basedir + "CylInitRe" + str(Re[i]);
        for j in range(0,times):
            ind = IC[i] + samp*j;
            numDig = len(str(ind));
            if numDig==3:
                FileName = FileNameBase + "/ibpm" + str(0) + str(0) + str(ind) + ".plt";
            elif numDig==4:
                FileName = FileNameBase + "/ibpm" + str(0) + str(ind) + ".plt";
            elif numDig==5:
                FileName = FileNameBase + "/ibpm" + str(ind) + ".plt";
                # Scan in data
                print "PROCESSING file" + str(iter+1) + " " + FileName + "...\n"
                data = IBPMData();
                data = readPltData(FileName);
                # Append data as column of POD matrix
                if ((i==0) & (j==0)):
                    X = np.zeros((data.nx*data.ny,snapshots));
                X[0:data.nx*data.ny,iter] = data.Z;
                iter = iter+1;
    # Output POD data matrix to file
    np.savetxt('PODMatrix.dat',X);
elif LoadRead == "READ":
    print "READING DATA MATRIX FROM FILE..."
    X = np.loadtxt('PODMatrix.dat');
    ind = 58000;
    FileName = basedir + "CylInitRe" + str(Re[0]) + "/ibpm" + str(ind) + ".plt";
    data = IBPMData();
    data = readPltData(FileName);
nx = data.nx; ny = data.ny;
sizeX = np.size(X,0);

# **********************************************************
# Global POD, UQ on DMD eigenvalues/vectors
# **********************************************************

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
    LAM[i] = (1.0/M)*mu[i];
    contourPlot(data.X,data.Y,PHI[0:sizeX,i]);
# Project data onto POD modes
coeff = np.zeros((M,snapshots));
for i in range(0,snapshots):
    for j in range(0,M):
        coeff[j,i] = np.dot(np.transpose(PHI[0:sizeX,j]),X[0:sizeX,i]);
for i in range(0,runs):
    indStart = times*i;
    indEnd = times*(i+1);
    plt.figure();
    for j in range(0,M):
        plt.plot(np.linspace(1,times,times),coeff[j,indStart:indEnd],marker="o");
# DMD
coeffDMD = np.zeros((M,snapshots/M));
eigDMD = np.zeros((M,runs),dtype=complex);
vecDMD = np.zeros((M,M,runs),dtype=complex);
for i in range(0,runs):
    indStart=times*i;
    indEnd=times*(i+1);
    A = np.dot(coeff[0:M,indStart+1:indEnd],np.linalg.pinv(coeff[0:M,indStart:indEnd-1]));
    coeffDMD[0:M,0] = coeff[0:M,indStart];
    coeffDMD[0:M,1:snapshots/M] = np.dot(A,coeff[0:M,indStart:indEnd-1]);
    muDMD,vDMD = np.linalg.eig(A);
    eigDMD[:,i] = muDMD;
    vecDMD[:,:,i] = vDMD;
    plt.figure();
    for j in range(0,M):
        plt.plot(np.linspace(1,times,times),coeff[j,indStart:indEnd],marker="o",color='blue');
        plt.plot(np.linspace(1,times,times),coeffDMD[j,0:snapshots/M],marker="x",color='red');
eigC = np.log(eigDMD)/DT;
# Separate eigenvalues into groups based on frequency
omegaHIGH = 1.5;
omegaLOW = 0.5;
eigHIGH = np.zeros(M);
eigLOW = np.zeros(M);
iterHIGH = 0;
iterLOW = 0;
eigC = np.reshape(eigC,M*runs,1);
for i in range(0,M*runs):
    if np.imag(eigC[i]) > omegaHIGH:
        eigHIGH[iterHIGH] = np.imag(eigC[i]);
        iterHIGH += 1;
    elif ((np.imag(eigC[i]) > omegaLOW) & (np.imag(eigC[i]) < omegaHIGH)):
        eigLOW[iterLOW] = np.imag(eigC[i]);
        iterLOW += 1;
# PCE on DMD eigenvalue distribution
weights = [0.1185, 0.2393, 0.2844, 0.2393, 0.1185];
xiQ = np.zeros(runs);
coeffHIGH = np.zeros(runs);
coeffLOW = np.zeros(runs);
for i in range(0,runs):
    xiQ[i] = (2.0/(endPts[np.size(endPts)-1]-endPts[0]))*(Re[i]-endPts[0]) - 1;
coeffHIGH = LegendreProjection(eigHIGH,xiQ,weights);
coeffLOW = LegendreProjection(eigLOW,xiQ,weights);
# PCE on DMD eigenvector distribution


# Sample, plot histogram of statistics
nsamps = 50000;
samps = np.random.uniform(-1,1,nsamps);
distHIGH = np.zeros(nsamps);
distLOW = np.zeros(nsamps);
for i in range(0,nsamps):
    for j in range(0,runs):
        distHIGH[i] += coeffHIGH[j]*Legendre(samps[i],j);
        distLOW[i] += coeffLOW[j]*Legendre(samps[i],j);
n, bins, patches = plt.hist(distHIGH, 25, normed=1, facecolor='b', alpha=0.75);
figure();
n2, bins2, patches2 = plt.hist(distLOW, 25, normed=1, facecolor='b', alpha=0.75);

# **********************************************************
# POD separately for each simulation; do UQ on POD modes
# **********************************************************

# Extract POD modes for each Re number
print "COMPUTING POD...\n";
M = 5;
MEANS = np.zeros((sizeX,runs));
POD = np.zeros((sizeX,M,runs));
for i in range(0,runs):
    ind1 = times*i;
    ind2 = times*(i+1);
    Xtmp = X[:,ind1:ind2];
    Xmean = np.sum(Xtmp,axis=1)/times;
    Xtmp = Xtmp - np.transpose(np.tile(Xmean,[times,1]));
    XTX = np.dot(np.transpose(Xtmp),Xtmp);
    mu,v = np.linalg.eig(XTX);
    MEANS[:,i] = Xmean;
    for j in range(0,M):
        POD[:,j,i] = 1/np.sqrt(mu[j])*np.dot(Xtmp,v[:,j]);
# Do PCE on the actual POD modes
weights = [0.1185, 0.2393, 0.2844, 0.2393, 0.1185];
xiQ = np.zeros(runs);
for i in range(0,runs):
    xiQ[i] = (2.0/(endPts[np.size(endPts)-1]-endPts[0]))*(Re[i]-endPts[0]) - 1;
Q = runs;
LEGcoeffMEAN = np.zeros((sizeX,Q));
LEGcoeffPOD = np.zeros((sizeX,Q,M));
for j in range(0,sizeX):
    LEGcoeffMEAN[j,:] = LegendreProjection(MEANS[j,:],xiQ,weights);
for j in range(0,Q):
    for k in range(0,sizeX):
        LEGcoeffPOD[k,j,:] = LegendreProjection(POD[k,j,:],xiQ,weights);
# Compute variance
VarianceMEAN = np.zeros((sizeX));
VariancePOD = np.zeros((sizeX,M));
for i in range(0,sizeX):
    VarianceMEAN[i] += np.sum(np.power(LEGcoeffMEAN[i,1:Q],2));
for i in range(0,sizeX):
    for j in range(0,M):
        VariancePOD[i,j] += np.sum(np.power(LEGcoeffPOD[i,j,1:Q],2));
# Compute reconstructions
REC = np.zeros((sizeX));
mode = 0;
xi = xiQ[4];
for i in range(0,sizeX):
    for j in range(0,Q):
        REC[i] += LEGcoeffPOD[i,mode,j]*Legendre(xi,j);

# Wait for user to manually end program
_ = raw_input("Press [enter] to continue.");
