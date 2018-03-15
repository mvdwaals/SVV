import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os
import codecs
from datastationary import *
from dataconst import *
from dataweight import *
from funcv import *
from functhrust import *

n_test = len(data_not_si)
data_not_si_T = data_not_si.T

hp_ft = data_not_si_T[0]
Vc_kts = data_not_si_T[1]
Tmta_c = data_not_si_T[2]
FFl_lbhr = data_not_si_T[3]
FFr_lbhr = data_not_si_T[4]
mfu_lb = data_not_si_T[5]
alpha_deg = data_not_si_T[6]
alpha_rad = np.radians(alpha_deg)

empty_weight = empty_weight_lb * lb_kg
fuel_weight = fuel_weight_lb * lb_kg
m_tot = sum(person_weight) + empty_weight + fuel_weight

hp = hp_ft * ft_m
Vc = Vc_kts * kts_ms
m = m_tot - mfu_lb * lb_kg
Tmta = Tmta_c + c_k
FFl = FFl_lbhr * lbhr_kgs
FFr = FFr_lbhr * lbhr_kgs

W = m * g

# Reductions to Ve
p = fp(hp)
M = fM(p, Vc)
T = fT(M, Tmta)
a = fa(T)
Vt = fVt(M, a)
rho = frho(p, T)

# Lift coefficient
CL = 2 * W / (rho * Vt**2 * S)              # Lift coefficient [ ]

#Plotting and finding CLa
CLa, intercept, r_value, uu_p_value, uu_std_err = stats.linregress(alpha_rad,CL) # Lots of unused (uu_) values

linregress_x = np.array([min(alpha_rad), max(alpha_rad)])
linregress_y = intercept + CLa * linregress_x

CLalabel = 'CLa = '+str(round(CLa,3))+' [5.084], r = '+str(round(r_value,3))
print(CLalabel)

plt.plot(linregress_x, linregress_y, label = CLalabel)
plt.scatter(alpha_rad, CL)
plt.ylabel('CL [-]')
plt.xlabel('alpha [rad]')
plt.legend()
plt.savefig('graphclalpha.png')
plt.cla()
plt.clf()

Ttotal = fTtotal(T,n_test, hp, M, FFl, FFr)

CD = 2 * Ttotal / (rho * Vt**2 * S)
 
CL_sq = CL**2

#Plotting CLsq-CD, find e and CD0
slope, CD0, r_value, uu_p_value, uu_std_err = stats.linregress(CL_sq,CD) # Lots of unused (uu_) values

linregress_x = np.array([min(CL_sq), max(CL_sq)])
linregress_y = CD0 + slope * linregress_x

oswald = 1 / (pi * A * slope)

CClabel = 'CD0 = '+str(round(CD0,3))+' [0.04], e = '+str(round(oswald,3))+'[0.8], r = '+str(round(r_value,3))

print(CClabel)

plt.plot(linregress_x, linregress_y, label=CClabel)
plt.scatter(CL_sq, CD)
plt.ylabel('CD [-]')
plt.xlabel('CL^2 [-]')
plt.legend()
plt.savefig('graphcl2cd.png')
plt.cla()
plt.clf()

#Plotting CL-CD
plt.scatter(CL, CD)
plt.ylabel('CD [-]')
plt.xlabel('CL [-]')
plt.savefig('graphclcd.png')
plt.cla()
plt.clf()

#Plotting CD-a
CDa, intercept, r_value, uu_p_value, uu_std_err = stats.linregress(alpha_rad,CD) # Lots of unused (uu_) values

linregress_x = np.array([min(alpha_rad), max(alpha_rad)])
linregress_y = intercept + CDa * linregress_x

CDalabel = 'CDa = '+str(round(CDa,3))+', r = '+str(round(r_value,3))
print(CDalabel)

plt.plot(linregress_x, linregress_y, label = CDalabel)
plt.scatter(alpha_rad, CD)
plt.ylabel('CD [-]')
plt.xlabel('alpha [rad]')
plt.legend()
plt.savefig('graphcdalpha.png')
plt.cla()
plt.clf()

print('Graphs exported.')

