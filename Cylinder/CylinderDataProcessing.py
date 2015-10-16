import sys
import numpy as np
import pylab
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

# Read file data
FileName = "Cylinder/CylInitRe97.6545/ibpm10000.bin";
data = IBPMData();
data = Bin2Dat(FileName);
contourPlot(data.X,data.Y,data.Z);
