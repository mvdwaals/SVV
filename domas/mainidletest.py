import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os
import codecs
from datatrim import *
from dataconst import *
from dataweight import *
from funcv import *

n_tests = len(data_not_si)
data_not_si_T = data_not_si.T

hp_ft = data_not_si_T[0]
Vc_kts = data_not_si_T[1]
mfu_lb = data_not_si_T[8]
Tmta_c = data_not_si_T[9]

empty_weight = empty_weight_lb * lb_kg
fuel_weight = fuel_weight_lb * lb_kg
m_tot = sum(person_weight) + empty_weight + fuel_weight

w_nose = 1000  #nose payload weight, assume smt
w_aftcabin =  1000 #aft cabin payload weight, assume smt

m = m_tot - mfu_lb * lb_kg

#parameters
cg_bem = 300.21   #cg basic empty mass in inches, CHANGE
MAC = 80.98       #mean aerodynamic chord
X_LEMAC = 261.56  #coordinate leading edge mean aerodynamic chord
xcg_int = [131,131,214,214,251,251,288,288,170,74,321]
xcg_int_moved = [131,131,131,214,214,251,251,288,288,170,74,321]

results = []
times = []
fuelusedleft = []
fuelusedright = []
fuelused = []

f= open('time.txt',"r")
for line in f:
    results.append(line.split())

for i in range(len(results)):
    for j in range(len(results[i])):
        times.append(results[i][j])
results = []

g= open('Fuelusedleft.txt',"r")
for line in g:
    results.append(line.split())

for i in range(len(results)):
    #for j in range(0,len(results[i])):
    fuelusedleft.append(results[i][0])
results = []

h= open('Fuelusedright.txt',"r")
for line in h:
    results.append(line.split())

for i in range(len(results)):
    #for j in range(0,len(results[i])):
    fuelusedright.append(results[i][0])

for i in range(len(fuelusedleft)):
    fuelused.append(float(fuelusedleft[i]) + float(fuelusedright[i]))


person_weight.extend([w_nose* lb_kg,w_aftcabin* lb_kg])
person_weight_moved.extend([w_nose* lb_kg,w_aftcabin* lb_kg])

w_int = person_weight
w_int_moved = person_weight_moved

xcg = cg_bem*in_m
w_tot = empty_weight + fuel_weight + sum(person_weight)  #total weight
w = w_tot   #set weight variable
xbarcg = []   #set list for getting cg positions
weight = []   #set list for getting weights
#xbarcg.append((xcg - X_LEMAC*in_m)/(MAC*in_m))  #set initial data point 
weight.append(w)

for i in range(len(times)): 

    if float(times[i]) < 3064 and float(times[i]) > 3162:
        for i in range(len(xcg_int)):
            mom_int = xcg_int[i]*in_m*w_int[i]*lb_kg
            w_tot_int = w_int[i]*lb_kg
    else:
        for i in range(len(xcg_int_moved)):
            mom_int = xcg_int_moved[i]*in_m*w_int_moved[i]*lb_kg
            w_tot_int = w_int_moved[i]*lb_kg

    xbarcg.append(((xcg*w + mom_int)/(w + w_tot_int) - X_LEMAC*in_m)*100/(MAC*in_m))
    xcg = ((xcg*w + mom_int)/(w + w_tot_int) - X_LEMAC*in_m)*100/(MAC*in_m)
    w = w_tot - fuelused[i]
    weight.append(w)

plt.plot(times,xbarcg)
plt.show()
print(fuelusedleft[4])
print(fuelusedright[4])
print(times[3162])
