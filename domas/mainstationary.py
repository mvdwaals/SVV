import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os
import codecs
from datastationary import *
from dataconst import *
from dataweight import *
from funcv import *

n_tests = len(data_not_si)
data_not_si_T = data_not_si.T

hp_ft = data_not_si_T[0]
Vc_kts = data_not_si_T[1]
Tmta_c = data_not_si_T[2]
FFl_lbhr = data_not_si_T[3]
FFr_lbhr = data_not_si_T[4]
mfu_lb = data_not_si_T[5]
alpha_deg = data_not_si_T[6]

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
CLa, intercept, uu_r_value, uu_p_value, uu_std_err = stats.linregress(alpha_deg,CL) # Lots of unused (uu_) values

linregress_x = np.array([min(alpha_deg), max(alpha_deg)])
linregress_y = intercept + CLa * linregress_x

CLalabel = 'CLa = '+str(CLa)
print(CLalabel)

plt.plot(linregress_x, linregress_y, label = CLalabel)
plt.scatter(alpha_deg, CL)
plt.ylabel('CL [-]')
plt.xlabel('alpha [deg]')
plt.legend()
plt.savefig('graphclalpha.png')
plt.cla()
plt.clf()

# CD calculations
T_ISA = T0 + labda * hp
DeltaT = T - T_ISA
with codecs.open('matlab.dat','w',encoding='utf8') as f:
    for i in range(n_tests):
        str1 = str(hp[i])+'\t'
        str2 = str(M[i])+'\t'
        str3 = str(DeltaT[i])+'\t'
        str4 = str(FFl[i])+'\t'
        str5 = str(FFr[i])+'\r\n'
        strtot = str1+str2+str3+str4+str5
        f.write(strtot)
    f.close()

os.system('java -jar thrust.jar')

#dummy = input("Press any key when thrust.dat is updated.")

infile = np.genfromtxt('thrust.dat').T
T1 = infile[0]
T2 = infile[1]
Ttotal = T1 + T2

CD = 2 * Ttotal / (rho * Vt**2 * S)
 
CL_sq = CL**2

#Plotting CLsq-CD, find e and CD0
slope, CD0, uu_r_value, uu_p_value, uu_std_err = stats.linregress(CL_sq,CD) # Lots of unused (uu_) values

linregress_x = np.array([min(CL_sq), max(CL_sq)])
linregress_y = CD0 + slope * linregress_x

oswald = 1 / (pi * A * slope)

CClabel = 'CD0 = '+str(CD0)+', e = '+str(oswald)

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
#p = np.polyfit(CD,CL,2)
#x = np.linspace(min(CL),max(CL),50)
#y = p[0]*x**2 + p[1]*x + p[2]

#plt.plot(x,y)
plt.plot(CD, CL, 's-')
plt.ylabel('CD [-]')
plt.xlabel('CL [-]')
plt.savefig('graphclcd.png')
plt.cla()
plt.clf()


#Plotting CD-a
CDa, intercept, uu_r_value, uu_p_value, uu_std_err = stats.linregress(alpha_deg,CD) # Lots of unused (uu_) values

linregress_x = np.array([min(alpha_deg), max(alpha_deg)])
linregress_y = intercept + CDa * linregress_x

CDalabel = 'CDa = '+str(CDa)
print(CDalabel)

plt.plot(linregress_x, linregress_y, label = CDalabel)
plt.scatter(alpha_deg, CD)
plt.ylabel('CD [-]')
plt.xlabel('alpha [deg]')
plt.legend()
plt.savefig('graphcdalpha.png')
plt.cla()
plt.clf()

print('Graphs exported.')

