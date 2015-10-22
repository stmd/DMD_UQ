import sys
import numpy as np
import pylab
from __future__ import division
# Custom auxiliary functions
import Legendre as Legendre

# **************************************************
# Discrete projection onto Legendre polynomials
# **************************************************
def LegendreProjection(fQ,xiQ,wQ):
    Q = np.size(wQ);
    PROJ = 0;
    for i in range(0,Q):
        LEG = Legendre(xiQ,i);
        NORM = 2/(2*i+1);
        PROJ = PROJ + LEG*fQ*wQ[i]/NORM;
    return PROJ;


# **************************************************
# Subroutine to return L_n(x)
# **************************************************
def Legendre(x,n):
    if n==0:
        return 1;
    elif n==1:
        return x;
    else:
        LEG = (2.0*n-1)/n*x*Legendre(x,n-1) - (n-1)/n*Legendre(x,n-2);
        return LEG;
