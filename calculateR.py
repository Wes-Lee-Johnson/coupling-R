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

