import sys
import numpy as np
import pylab

# **********************************************************
# Struct for handling data IO
# **********************************************************
class IBPMData:
    def __init__(self):
        self.nx = 0;
        self.ny = 0;
        self.ngrid = 0;
        self.dx = 0;
        self.x0 = 0;
        self.y0 = 0;
        self.numPoints = 0;
        self.q = 0;
        self.FX = 0;
        self.FY = 0;
        self.timestep = 0;
        self.time = 0;
        self.X = 0;
        self.Y = 0;
        self.U = 0;
        self.V = 0;
        self.Z = 0;
