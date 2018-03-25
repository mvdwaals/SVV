import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os
import codecs
from datatrim import *
from dataconst import *
from dataweight import *
from funcv import *
from funcfuelarm import findarm
from funcfileimport import times, fuelusedleft, fuelusedright

#measured times and fuel used
t = times()
ful = fuelusedleft()
fur = fuelusedright()
fut = []

for i in range(len(ful)):
    fut.append(ful[i]+fur[i])

#unit conversion
n_tests = len(data_not_si)
data_not_si_T = data_not_si.T

hp_ft = data_not_si_T[0]
Vc_kts = data_not_si_T[1]
mfu_lb = data_not_si_T[8]
Tmta_c = data_not_si_T[9]

#weights
w_nose = 0  #nose payload weight
w_aftcabin =  220*lb_kg #aft cabin payload weight
person_weight.extend([w_nose,w_aftcabin])
person_weight_moved.extend([w_nose,w_aftcabin])
w_pl = person_weight  #payload weight
w_pl_moved = person_weight_moved #payload weight when Myk moved
w_oe = empty_weight_lb * lb_kg  #operational empty weight
w_f = fuel_weight_lb * lb_kg #fuel weight
w_zpl = w_oe + w_f      #zero payload weight
w_tot = sum(w_pl) + w_oe + w_f   #total weight

#parameters
cg_bem = 292.18*in_m      #cg basic empty mass 
MAC = 80.98*in_m       #mean aerodynamic chord
X_LEMAC = 261.45*in_m  #coordinate leading edge mean aerodynamic chord
coordinates = [131,131,214,214,251,251,288,288,170,74,321]      
xcg_pl = []       #payload positions
coordinates_moved = [131,131,131,214,214,251,251,288,288,170,74,321]    
xcg_pl_moved = []  #payload positions when Myk moved

for i in range(len(coordinates)):
    xcg_pl.append(coordinates[i]*in_m)

for i in range(len(coordinates_moved)):
    xcg_pl_moved.append(coordinates_moved[i]*in_m)

#loop
xcg = cg_bem
xbarcg = []   #set list for getting normalized cg positions (with MAC)
xcgtot = []   #set list for getting absolute cg positions
index = []    #find index movement
t_init_myke = 3064   #time in seconds that Myke moved forward
t_end_myke = 3162    #time in seconds that Myke moved aft
    
for i in range(len(t)):    
    mom_pl = 0
    w_tot_pl = 0
    if t[i] < t_init_myke or t[i] > t_end_myke:
        for j in range(len(xcg_pl)):
            mom_pl += xcg_pl[j]*w_pl[j]
            w_tot_pl += w_pl[j]
    else:
        index.append(i)
        for j in range(len(xcg_pl_moved)):
            mom_pl += xcg_pl_moved[j]*w_pl_moved[j]
            w_tot_pl += w_pl_moved[j]          
    xbarcg.append(((xcg*w_oe + mom_pl + (w_f - fut[i])*(findarm(w_f - fut[i])+xcg))/(w_oe + w_tot_pl + (w_f - fut[i]))-X_LEMAC)*100/MAC)
    xcgtot.append((xcg*w_oe + mom_pl + (w_f - fut[i])*(findarm(w_f - fut[i])+xcg))/(w_oe + w_tot_pl + (w_f - fut[i])))

#find delta xcg
def fDeltaxcg():
    
    leftjump = xcgtot[index[0]] - xcgtot[index[0] - 1]
    rightjump = xcgtot[index[len(index)-1]+1] - xcgtot[index[len(index)-1]]

    return (abs(leftjump) + abs(rightjump))/2, (abs(leftjump) + abs(rightjump))/2*100/MAC

#plotting
plt.xlabel("time [s]")
plt.ylabel("x_cg [%MAC]")
plt.plot(t,xbarcg)
plt.grid()
plt.savefig('graphcgt.png')

