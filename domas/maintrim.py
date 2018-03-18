import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os
import codecs
from datatrim import *
from dataconst import *
from dataweight import *
from funcv import *
from funccg import *
from functhrust import *

def avg(lst): return sum(lst)/len(lst)

n_r = 3

data_not_si = myk_not_si
n_tests = len(data_not_si)
data_not_si_T = data_not_si.T

hp_ft = data_not_si_T[0]
Vc_kts = data_not_si_T[1]
alpha_deg = data_not_si_T[2]
deltae_deg = data_not_si_T[3]
Fe = data_not_si_T[5]
FFl_lbhr = data_not_si_T[6]
FFr_lbhr = data_not_si_T[7]
mfu_lb = data_not_si_T[8]
Tmta_c = data_not_si_T[9]

empty_weight = empty_weight_lb * lb_kg
fuel_weight = fuel_weight_lb * lb_kg
m_tot = sum(person_weight) + empty_weight + fuel_weight

hp = hp_ft * ft_m
Vc = Vc_kts * kts_ms
alpha_rad = np.radians(alpha_deg)
deltae_rad = np.radians(deltae_deg)
m = m_tot - mfu_lb * lb_kg
FFl = FFl_lbhr * lbhr_kgs
FFr = FFr_lbhr * lbhr_kgs
Tmta = Tmta_c + c_k

W = m * g

p = fp(hp) #hp
M = fM(p, Vc) #Vc
T = fT(M, Tmta) #Tmta
a = fa(T)
Vt = fVt(M, a)
rho = frho(p, T)
Ve = fVe(Vt, rho)
Vetilde = (fVetilde(Ve, Ws, W))[:-1] #Ignore the last measurement (after deltacg)

CN = 2 * W / (rho * Vt**2 * S)
CN = CN[-2::] #Assumption!!!

Deltaxcg = fDeltaxcg()[1] 
ddeltae_dalpha, uu_intercept, r_value, uu_p_value, uu_std_err = stats.linregress(alpha_rad,deltae_rad) # Lots of unused (uu_) values
Cmdelta = -CN * Deltaxcg / Deltadeltae_rad / c
#print(Cmdelta)
Cmdelta = avg(Cmdelta)

Cmalpha = -ddeltae_dalpha * Cmdelta
#print(Cmalpha)

print('Cmdelta = '+str(round(Cmdelta,n_r))+', Cmalpha = '+str(round(Cmalpha,n_r))+' with r = '+str(round(r_value,n_r)))

Ttotal = fTtotal(T,n_tests, hp, M, FFl, FFr)
Ts = fTtotal([T0], 1, [0], [0], [mdotfs], [mdotfs]) #Needs checking
Tcs = Ts * 2 / (rho * Vt**2 * S)
Tc = Ttotal * 2 / (rho * Vt**2 * S)
deltastareeq = (deltae_rad - CmTc * (Tcs - Tc) / Cmdelta)[:-1] #Ignore the last measurement (after deltacg)
Fstareaer = (Fe * Ws / W)[:-1] #Ignore the last measurement (after deltacg)

arr = np.array([Vetilde, deltastareeq, Fstareaer]).T
arr = arr[arr[:,0].argsort()]
arr = arr.T

plt.plot(arr[0], arr[1], 's-')
plt.xlabel('Ṽ_e [m/s]')
plt.ylabel('delta*_eq [rad]')
plt.savefig('graphdeltastar.png')
plt.cla()
plt.clf()

plt.plot(arr[0], arr[2], 's-')
plt.xlabel('Ṽ_e [m/s]')
plt.ylabel('F*_eq [N]')
plt.savefig('graphfstar.png')
plt.cla()
plt.clf()
#print(Vetilde)
#print(deltastareeq)
#print(Fstareaer)
