import sys
import numpy as np
import pylab
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter, ScalarFormatter, MaxNLocator
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure, axes, plot, xlabel, ylabel, title, grid, savefig, show

# **********************************************************
# Function to wrap contour plot details
# **********************************************************

def contourPlot(X,Y,Z):
    # Set figure parameters
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
    pylab.rcParams.update(params);
    plt.ion();
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
    minC = 0.5*Z.min();
    maxC = 0.5*Z.max();
    plt.scatter(X,Y,c=Z,cmap=cm.coolwarm,s=5,lw=0,vmin=minC,vmax=maxC);
    ax2 = plt.gca();
    ax2.set_aspect('equal')
    cbar = plt.colorbar(CS3, ticks=[0, 100])
    cbar.ax.set_yticklabels(["{0:.2f}".format(minC), "{0:.2f}".format(maxC)])

    ax.set_axisbelow(True)
    plt.gcf().subplots_adjust(left=0.16)
    plt.gcf().subplots_adjust(bottom=0.15)
    ax.xaxis.labelpad = 20
    plt.gcf().tight_layout()
    #plt.grid(None)
    pylab.savefig('contourPlot.eps',bbox_inches=0)
    #plt.gca().tight_layout()
    plt.draw()
    plt.show()
