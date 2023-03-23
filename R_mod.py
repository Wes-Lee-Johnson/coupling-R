"""
Author: Wes Johnson
Date:   03112
Purpose:
        This Script analyzes the simulation by loading the pickled 
        mode analyist code and saved '.npy' files containing the 
        trajectories of the ions in 6D phase space. 
        Contains functions to calculate the energy in rotating frame 
        indepedent of mode analysis.

        Analyzes data for 4ion parameter scan over offset and width 
        of perp. laser.
"""
#cc05ccc10ccc15ccc20ccc25ccc30ccc35ccc40ccc45ccc50ccc55ccc60ccc65ccc70cc
# imports 
import numpy as np
import matplotlib.pyplot as plt 
from scipy import constants 

#cc05ccc10ccc15ccc20ccc25ccc30ccc35ccc40ccc45ccc50ccc55ccc60ccc65ccc70cc
#functions 

def calc_wz_b(wr,wc,b=0.0348): 
    return np.sqrt( (wr*wc - wr**2)/(b + 1/2) )

def calc_wz_R(wr,wc,R=5): 
    a   = (wc - 2*wr)*(R + 1)/(R - 1)
    return np.sqrt( (1/2)*(wc**2 - a**2) )


def findMinDiff(x,y,z):
    diff    = np.abs(x-y)
    notnan  = np.logical_not((np.isnan(diff)))
    mindex  = diff==np.min(diff[notnan])
    return  1/2*(x[mindex]+y[mindex]),z[mindex]

def calc_wr(wc,R,b): 
    RpRm= ((R + 1)/(R - 1))**2  
    C   = 1/2*wc**2*(1 - RpRm)
    B   = wc*(1/(b + 1/2) - 2*RpRm)
    A   = B/wc 
    return (B + np.sqrt(B**2 - 4*A*C))/(2*A)

def calc_wpm(wz,wc):
    dw  = np.sqrt(wc**2 - 2*wz**2)
    wp  = 1/2*(wc + dw)
    wm  = 1/2*(wc - dw)
    return wp,wm 

def calc_R(wz,wc,wr): 
    wp,wm   = calc_wpm(wz,wc) 
    return abs((wp - wr)/(wm - wr))

def calc_b(wz,wc,wr): 
    beta = (wr*wc - wr**2)/(wz**2) - 1/2
    return beta 

