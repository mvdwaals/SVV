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

n_r = 3 # Round to number of digits
usealldata = 2 # 0 = manual data, 1 = manual data + trim data, 2 = matlab data

if usealldata == 2: # Read all data from dataflat.txt (Configure this file using maindata.py)
    print("Using matlab data.")
    data_not_si = np.genfromtxt('dataflat.txt')
    data_not_si_T = data_not_si
    hp_ft = data_not_si_T[3]
    Vc_kts = data_not_si_T[6]
    Tmta_c = data_not_si_T[5]
    FFl_lbhr = data_not_si_T[1]
    FFr_lbhr = data_not_si_T[2]
    mfu_lb = data_not_si_T[4]
    alpha_deg = data_not_si_T[0]
    n_test = len(hp_ft)
    saved = True

else: # Read data from manually recorded data
    print("/!\ Not using matlab data.")
    if usealldata == 1: data_not_si = np.concatenate((data_not_si, trim_not_si))
    saved = False
    data_not_si_T = data_not_si.T
    n_test = len(data_not_si)
    hp_ft = data_not_si_T[0]
    Vc_kts = data_not_si_T[1]
    Tmta_c = data_not_si_T[2]
    FFl_lbhr = data_not_si_T[3]
    FFr_lbhr = data_not_si_T[4]
    mfu_lb = data_not_si_T[5]
    alpha_deg = data_not_si_T[6]

print('Loading manual data.')
datasaved = np.genfromtxt('datamanual.txt')
alphas = datasaved[0]
CLs = datasaved[1]
CDs = datasaved[2]


# Convert all to SI units
empty_weight = empty_weight_lb * lb_kg
fuel_weight = fuel_weight_lb * lb_kg
m_tot = sum(person_weight_value) + empty_weight + fuel_weight

alpha_rad = np.radians(alpha_deg)
hp = hp_ft * ft_m
Vc = Vc_kts * kts_ms
m = m_tot - mfu_lb * lb_kg
Tmta = Tmta_c + c_k
FFl = FFl_lbhr * lbhr_kgs
FFr = FFr_lbhr * lbhr_kgs

W = m * g

alpha = alpha_rad # Choose between alpha_rad and alpha_deg

# Intermediate steps in reductions to Ve, Ve itself is not used
p = fp(hp)
M = fM(p, Vc)
T = fT(M, Tmta)
a = fa(T)
Vt = fVt(M, a)
rho = frho(p, T)

#Re and M range
mu_air = labda_air * T**(3/2) / (T + C_air)
Re = rho * Vt * c / mu_air

Mrange = str(round(min(M),n_r))+' - '+str(round(max(M),n_r))
Rerange = str(int(min(Re)))+' - '+str(int(max(Re)))

# Lift coefficient
CL = 2 * W / (rho * Vt**2 * S)              # Lift coefficient [ ]

#Plotting and finding CLa by linear regression
CLa, intercept, r_value, uu_p_value, uu_std_err = stats.linregress(alpha,CL) # Lots of unused (uu_) values

linregress_x = np.array([min(alpha), max(alpha)])
linregress_y = intercept + CLa * linregress_x

CLalabel = '$C_{L_a}$ = '+str(round(CLa,n_r))+' [5.084], $r^2$ = '+str(round(r_value**2,n_r))
print(CLalabel)
plt.plot(linregress_x, linregress_y, label = CLalabel)
plt.scatter(alpha, CL, label = 'Automatically recorded')
plt.scatter(alphas, CLs, label = 'Manually recorded')
plt.title('$C_L / \alpha$ at clean configuration,\n Mach range = '+Mrange+', Re range = '+Rerange)
plt.ylabel('$C_L$ [-]')
plt.xlabel('\alpha [rad]')
plt.grid()
plt.legend()
plt.savefig('graphclalpha.png')
plt.cla()
plt.clf()

# Calculate thrust using provided Java exectable
print('Running Java program.')
Ttotal = fTtotal(T,n_test, hp, M, FFl, FFr)
print('Java program finished.')

CD = 2 * Ttotal / (rho * Vt**2 * S)
 
CL_sq = CL**2

#Plotting CLsq-CD, and find e and CD0, using linear regression
slope, CD0, r_value, uu_p_value, uu_std_err = stats.linregress(CL_sq,CD) # Lots of unused (uu_) values

linregress_x = np.array([min(CL_sq), max(CL_sq)])
linregress_y = CD0 + slope * linregress_x

oswald = 1 / (pi * A * slope)

CClabel = '$C_{D_0}$ = '+str(round(CD0,n_r))+' [0.04], $e$ = '+str(round(oswald,n_r))+' [0.8], r^2 = '+str(round(r_value**2,n_r))

print(CClabel)

plt.plot(linregress_x, linregress_y, label=CClabel)
plt.scatter(CL_sq, CD, label = 'Automatically recorded')
plt.scatter(CLs**2, CDs, label = 'Manually recorded')
plt.title('$C_D / C_L^2$ at clean configuration,\n Mach range = '+Mrange+', Re range = '+Rerange)
plt.ylabel('$C_D$ [-]')
plt.xlabel('$C_L^2$ [-]')
plt.legend()
plt.grid()
plt.savefig('graphcl2cd.png')
plt.cla()
plt.clf()

#Plotting CL-CD
fit = np.polyfit(CL,CD,2)
x = np.linspace(min(CL),max(CL))
plt.plot(x,fit[0]*x**2 + fit[1] *x + fit[2])
plt.scatter(CL, CD, label = 'Automatically recorded')
plt.scatter(CLs, CDs, label = 'Manually recorded')
plt.title('$C_D / C_L$ at clean configuration,\n Mach range = '+Mrange+', Re range = '+Rerange)
plt.ylabel('$C_D$ [-]')
plt.xlabel('$C_L$ [-]')
plt.legend()
plt.grid()
plt.savefig('graphclcd.png')
plt.cla()
plt.clf()

#Plotting CD-a
CDa, intercept, r_value, uu_p_value, uu_std_err = stats.linregress(alpha,CD) # Lots of unused (uu_) values

linregress_x = np.array([min(alpha), max(alpha)])
linregress_y = intercept + CDa * linregress_x

CDalabel = '$C_{D_a}$ = '+str(round(CDa,n_r))+', $r^2$ = '+str(round(r_value**2,n_r))
print(CDalabel)

#plt.plot(linregress_x, linregress_y, label = CDalabel)
#plt.legend()
plt.scatter(alpha, CD, label = 'Automatically recorded')
plt.scatter(alphas, CDs, label = 'Manually recorded')
plt.title('$C_D / \alpha$ at clean configuration,\n Mach range = '+Mrange+', Re range = '+Rerange)
plt.ylabel('$C_D$ [-]')
plt.xlabel('$\alpha$ [rad]')
plt.legend()
plt.grid()
plt.savefig('graphcdalpha.png')
plt.cla()
plt.clf()

if not saved:
    np.savetxt('datamanual.txt', np.array([alpha, CL, CD]))
    print('Manual data saved.')

print('Graphs exported.')

