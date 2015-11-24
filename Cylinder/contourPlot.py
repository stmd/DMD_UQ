import sys
import numpy as np
import pylab
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter, ScalarFormatter, MaxNLocator
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, axes, plot, xlabel, ylabel, title, grid, savefig, show
from mpl_toolkits.axes_grid1 import make_axes_locatable

# **********************************************************
# Function to wrap contour plot details
# **********************************************************

def contourPlot(X,Y,Z):
    fig = plt.figure()
    ax = plt.gca()
    # Colorbar code ************
    min, max = (0,100)
    step = 1
    # Using contourf to provide my colorbar info, then clearing the figure
    ZZ = [[0,0],[0,0]]
    Xtmp = np.array([np.min(X),np.max(X)]);
    Ytmp = np.array([np.min(Y),np.max(Y)]);
    XX,YY = np.meshgrid(Xtmp,Ytmp);
    levels = range(min,max+step,step)
    CS3 = plt.contourf(XX,YY,ZZ, levels, cmap=cm.coolwarm,aspect='equal')
    plt.gca().set_xlim([np.min(X),np.max(X)])
    plt.gca().set_ylim([np.min(Y),np.max(Y)])
    plt.clf()
    minC = 0.5*Z.min();
    maxC = 0.5*Z.max();
    maxC = 2.5e-4;
    minC = 0;
    #minC = 0;
    #maxC = 1e-4;
    plt.scatter(X,Y,c=Z,cmap=cm.coolwarm,s=5,lw=0,vmin=minC,vmax=maxC);
    axes().set_xlim([0.5,14]);
    axes().set_ylim([-2.0,2.0]);
    ax2 = plt.gca();
    ax2.set_aspect('equal')
    cbar = plt.colorbar(CS3, ticks=[0, 100],shrink=0.7,aspect=15)
    cbar.ax.set_yticklabels(["{:.0E}".format(minC), "{:.0E}".format(maxC)])
    ax.set_axisbelow(True)
    ax.xaxis.labelpad = 20
    #plt.gca().xaxis.set_major_locator(plt.NullLocator())
    #plt.gca().yaxis.set_major_locator(plt.NullLocator())
    plt.gcf().tight_layout()
    #plt.grid(None)
    pylab.savefig('contourPlot.eps',bbox_inches=0)
    #plt.gca().tight_layout()
    plt.draw()
    plt.show()
