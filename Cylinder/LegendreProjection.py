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
    PROJ = np.zeros(Q);
    for i in range(0,Q):
        for j in range(0,Q):
            LEG = Legendre(xiQ[j],i);
            NORM = 2.0/(2*i+1);
            PROJ[i] += LEG*fQ[j]*wQ[j]/NORM;
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
