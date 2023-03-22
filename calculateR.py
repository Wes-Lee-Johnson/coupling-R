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

# parameters 
B   = .9145
m   = 40.078*constants.u #calcium
q   = constants.e 
wz  = 2*np.pi * 116e3
wc  = q/m*B
wr  = 2*np.pi * 23.5e3 
#cc05ccc10ccc15ccc20ccc25ccc30ccc35ccc40ccc45ccc50ccc55ccc60ccc65ccc70cc

# functions 
#def calc_wc(q,m,B): 
#    """unmodified cyclotron frequency
#        q   :: charge 
#        m   :: mass 
#        B   :: magnetic feild
#    """
#    return q/m*B 

def calc_wz_b(wr,wc=q/m*B,b=0.0348): 
    return np.sqrt( (wr*wc - wr**2)/(b + 1/2) )

def calc_wz_R(wr,wc=q/m*B,R=5): 
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

#cc05ccc10ccc15ccc20ccc25ccco30ccc35ccc40ccc45ccc50ccc55ccc60ccc65ccc70cc
b   = calc_b(wz,wc,wr)
print('R=',calc_R(wz,wc,wr))
print('b=',b)
print('\n')
wr  = calc_wr(wc,5,b)
print('fr=',wr/2/np.pi/1e3)
print('fz=',calc_wz_b(wr,wc,b)/2/np.pi/1e3)
print('\n')
#cc05ccc10ccc15ccc20ccc25ccc30ccc35ccc40ccc45ccc50ccc55ccc60ccc65ccc70cc



#cc05ccc10ccc15ccc20ccc25ccc30ccc35ccc40ccc45ccc50ccc55ccc60ccc65ccc70cc
wz  = 2*np.pi * 1.58e6 
m   = 9.012*constants.u 
B   = 4.4588
wc  = q/m*B 
wr  = 2*np.pi * 180e3

print('R=',calc_R(wz,wc,wr))
print('b=',calc_b(wz,wc,wr))
print('fr=',calc_wr(wc,5,calc_b(wz,wc,wr))/2/np.pi/1e3)
exit()
#cc05ccc10ccc15ccc20ccc25ccc30ccc35ccc40ccc45ccc50ccc55ccc60ccc65ccc70cc




#cc05ccc10ccc15ccc20ccc25ccc30ccc35ccc40ccc45ccc50ccc55ccc60ccc65ccc70cc

# plotting 
unit1   = 2*np.pi*1e6
unit2   = 2*np.pi*1e3
wrmax   = q/m/2*B
wrmin   = 0
n   = 10**6
wr  = np.linspace(0e0,q/m*B/2,1000)
wzb = calc_wz_b(wr)
wzR = calc_wz_R(wr)
fig, ax  = plt.subplots(figsize=(8,8));
ax.plot(wr/unit2,wzb/unit1
        ,label=r'$\beta=$'+'%1.1e'%b
        )
ax.plot(wr/unit2,wzR/unit1
        ,label=r'$R=$'+'%1.0f'%R
        )
wzS,wrS = findMinDiff(wzb,wzR,wr)
#ax.hlines(y     = wzS/unit1
#        ,xmin   = wrmin/unit2
#        ,xmax   = wrmax/unit2
#        ,color  = 'k'
#        ,label=r'$f_z=$'+'%1.2f[MHz]'%(wzS/unit1)
#        )
wrA = calc_wr(q/m*B,R,b)
wzA = calc_wz_R(wrA)
print(wrA,wzA)
print((wrS-wrA)/wrA,(wzS-wzA)/wzA)
ax.scatter(wrS/unit2,wzS/unit1
        ,color  = 'k'
        ,label  = r'$f_z,f_{wall}=$'+'%1.2f[MHz],'%(wzA/unit1)+' %1.0f[kHz]'%(wrA/unit2)
        ,s      = 100
        ,marker = 'X'
        ,zorder = 10
        )
ax.set_ylabel('Axial Frequency MHz')
ax.set_xlabel('Wall  Frequency kHz')
ax.legend(loc='lower right')
plt.show()

