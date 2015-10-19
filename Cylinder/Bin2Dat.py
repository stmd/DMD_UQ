import sys
import numpy as np
from IBPMData import IBPMData

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
        # Set up (x,y) coordinates (corresponding to finest grid)
        xf = x0 + (nx-1)*dx;
        yf = y0 + (ny-1)*dx;
        x = np.linspace(x0,xf,nx);
        y = np.linspace(y0,yf,ny);
        Z = np.squeeze(omega[0,0:nx,0:ny]);
        Z = np.transpose(Z);
        XX,YY = np.meshgrid(x,y);
        XX = np.transpose(XX);
        YY = np.transpose(YY);
        X = np.zeros(nx*ny); Y = np.zeros(nx*ny);
        iter = 0;
        for i in range(0,ny):
            for j in range(0,nx):
                X[iter] = XX[j,i];
                Y[iter] = YY[j,i];
                iter = iter+1;
        # Return struct containing information
        data = IBPMData();
        data.nx = nx;
        data.ny = ny;
        data.ngrid = ngrid;
        data.dx = dx;
        data.x0 = x0;
        data.y0 = y0;
        data.numPoints = numPoints;
        data.q = q;
        data.FX = FX;
        data.FY = FY;
        data.timestep = timestep;
        data.time = time;
        data.X = X;
        data.Y = Y;
        data.Z = Z;
        
        return data;

    finally:
        fid.close()

