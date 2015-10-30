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
# Subroutine to generate a perturbed cylinder geometry for IBPM
# **********************************************************

def generatePerturbedCylinder(R,numPts,theta):
    dTH = 2*np.pi/numPts;
    th = np.zeros(numPts);
    for i in range(0,numPts):
        th[i] = i*dTH;
    xQ = theta;
    startPerturb = xQ*(np.pi/180);
    endPerturb = (xQ+20)*(np.pi/180);
    height = 0.2*R;
    gridStart = np.argmin(np.abs(th-startPerturb));
    gridEnd = np.argmin(np.abs(th-endPerturb));
    gridMiddle = int (np.round((gridStart+gridEnd)/2.0));
    Perturb = np.zeros(numPts);
    xMiddle = (R+height)*np.cos(th[gridMiddle]);
    yMiddle = (R+height)*np.sin(th[gridMiddle]);
    xStart = R*np.cos(th[gridStart]);
    yStart = R*np.sin(th[gridStart]);
    xEnd = R*np.cos(th[gridEnd]);
    yEnd = R*np.sin(th[gridEnd]);
    dx1 = (xMiddle-xStart)/(gridMiddle-gridStart);
    dy1 = (yMiddle-yStart)/(gridMiddle-gridStart);
    dx2 = (xEnd-xMiddle)/(gridEnd-gridMiddle);
    dy2 = (yEnd-yMiddle)/(gridEnd-gridMiddle);
    x = R*np.cos(th);
    y = R*np.sin(th);
    for i in range(gridStart,gridMiddle):
        x[i] = xStart + dx1*(i-gridStart);
        y[i] = yStart + dy1*(i-gridStart);
    for i in range(gridMiddle,gridEnd):
        x[i] = xMiddle + dx2*(i-gridMiddle);
        y[i] = yMiddle + dy2*(i-gridMiddle);
    xy = np.transpose([x,y]);
    return xy;

# **********************************************************
# Generate several perturbed cylinder geometries for IBPM
# **********************************************************

#xi = [82.1110, 90.3844, 102.5000, 114.6156, 122.8890];
xi = [84.2219, 100.7689, 125.0000, 149.2311, 165.7781];
Q = np.size(xi);
R = 0.5;
numPts = 160;
basedir = "/home/adegenna/ibpmcontrol/ibpmcontrol/output/PerturbedCylinder/"
for i in range(0,Q):
    # Generate cylinders
    xy = generatePerturbedCylinder(R,numPts,xi[i]);
    outFile = basedir + "grid" + str(xi[i]) + ".dat";
    fp = open(outFile,"w");
    fp.write("%d\n" %(numPts));
    for j in range(0,numPts):
        fp.write("%f %f\n" %(xy[j,0],xy[j,1]));
    fp.close();
    #Qsub the jobs
    header1 = "body PerturbedCylinder" + str(xi[i]);
    header2 = "raw " + outFile;
    header3 = "motion Stationary 0 0 0 0";
    header4 = "end"
    fileGeom = "PerturbedCylinder" + str(xi[i]) + ".dat";
    fGeom = open(fileGeom,"w");
    fGeom.write(header1+"\n");
    fGeom.write(header2+"\n");
    fGeom.write(header3+"\n");
    fGeom.write(header4+"\n");
    fGeom.close();
    outDir = basedir + "Grid" + str(xi[i]);
    geomFile = basedir + fileGeom;
    ic = outDir + "/ibpm08000.bin"
    os.system("qsub -v OUTDIR=" + outDir + ",GEOM=" + geomFile + " " + basedir + "pbsbase.dat");
    
