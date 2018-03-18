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

results = []

f= open('time.txt',"r")
for line in f:
    results.append(line.strip( ).split(','))

#for i in range(len()):

person_weight.extend([w_nose* lb_kg,w_aftcabin* lb_kg])

w_int = person_weight

xcg = cg_bem*in_m
w = empty_weight + fuel_weight + sum(person_weight)  #set weight variable
xbarcg = []   #set list for getting cg positions
weight = []   #set list for getting weights
xbarcg.append((cg_oew - X_LEMAC)/MAC)  #set initial data point 
weight.append(w)

#for i in range(100):   #change to list 

    #if t < VALUE and t < VALUE:
        #for i in range(len(xcg_int)):
            #mom_int = xcg_int[i]*in_m*w_int[i]*lb_kg
            #w_tot_int = w_int[i]*lb_kg
    #else:
        #for i in range(len(xcg_int)):
            #mom_int = xcg_int[i]*in_m*w_int[i]*lb_kg
            #w_tot_int = w_int[i]*lb_kg

    #xbarcg.append(((xcg*w + mom_int)/(w + w_tot_int) - X_LEMAC*in_m)*100/(MAC*in_m))
    #xcg = ((xcg*w + mom_int)/(w + w_tot_int) - X_LEMAC*in_m)*100/(MAC*in_m)
    #w -= LISTFUELFLOW[i]
    #weight.append(w)

print(results)
