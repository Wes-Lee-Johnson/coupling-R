"""
Author: Wes Johnson
Date:   032323
Purpose:
        This script compares the trapping frequencies need to achieve 
        particular values of beta and R
"""
#cc05ccc10ccc15ccc20ccc25ccc30ccc35ccc40ccc45ccc50ccc55ccc60ccc65ccc70cc
# imports 
import numpy as np
import matplotlib.pyplot as plt 
from scipy import constants 
import R_mod 

#cc05ccc10ccc15ccc20ccc25ccc30ccc35ccc40ccc45ccc50ccc55ccc60ccc65ccc70cc

# parameters 
B   = .9145
m   = 40.078*constants.u #calcium
q   = constants.e 
wz  = 2*np.pi * 116e3
wc  = q/m*B
wr  = 2*np.pi * 23.5e3 

#cc05ccc10ccc15ccc20ccc25ccco30ccc35ccc40ccc45ccc50ccc55ccc60ccc65ccc70cc
#parameters for 40Ca+
b   = R_mod.calc_b(wz,wc,wr)
print('R=',R_mod.calc_R(wz,wc,wr))
print('b=',b)
print('\n')
wr  = R_mod.calc_wr(wc,5,b)
print('fr=',wr/2/np.pi/1e3)
print('fz=',R_mod.calc_wz_b(wr,wc,b)/2/np.pi/1e3)
print('\n')
#cc05ccc10ccc15ccc20ccc25ccc30ccc35ccc40ccc45ccc50ccc55ccc60ccc65ccc70cc



#cc05ccc10ccc15ccc20ccc25ccc30ccc35ccc40ccc45ccc50ccc55ccc60ccc65ccc70cc
#parameters for 9Be+
wz  = 2*np.pi * 1.58e6 
m   = 9.012*constants.u 
B   = 4.4588
wc  = q/m*B 
wr  = 2*np.pi * 180e3

print('R=',R_mod.calc_R(wz,wc,wr))
print('b=',R_mod.calc_b(wz,wc,wr))
print('fr=',R_mod.calc_wr(wc,5,R_mod.calc_b(wz,wc,wr))/2/np.pi/1e3)
exit()
#cc05ccc10ccc15ccc20ccc25ccc30ccc35ccc40ccc45ccc50ccc55ccc60ccc65ccc70cc
