import numpy as np
import csv
# User-defined auxiliary functions/structures
from IBPMData import IBPMData

# *****************************************************
# Function to read IBPM data from .plt file
# *****************************************************

def readPltData(filename):
    with open(filename,'r') as f:
        next(f); next(f); next(f); next(f); next(f); next(f);
        reader = csv.reader(f,delimiter=' ');
        x=[];y=[];u=[];v=[];omega=[];tmp=[];
        for xi,yi,ui,vi,omi,tmpi in reader:
            x.append(xi); 
            y.append(yi); 
            u.append(ui); 
            v.append(vi); 
            omega.append(omi); 
            tmp.append(tmpi);
    
    data = IBPMData();
    data.X = x;
    data.Y = y;
    data.U = u;
    data.V = v;
    data.Z = omega;
    # Fix this!
    data.nx = 799;
    data.ny = 199;

    return data;
