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
alpha_deg = data_not_si_T[2]
deltae_deg = data_not_si_T[3]
Fe = data_not_si_T[6]
mfu_lb = data_not_si_T[8]
Tmta_c = data_not_si_T[9]

empty_weight = empty_weight_lb * lb_kg
fuel_weight = fuel_weight_lb * lb_kg
m_tot = sum(person_weight) + empty_weight + fuel_weight

hp = hp_ft * ft_m
Vc = Vc_kts * kts_ms
alpha_rad = np.radians(alpha_deg)
deltae_rad = np.radians(deltae_rad)
m = m_tot - mfu_lb * lb_kg
Tmta = Tmta_c + c_k

W = m * g

p = fp(hp) #hp
M = fM(p, Vc) #Vc
T = fT(M, Tmta) #Tmta
a = fa(T)
Vt = fVt(M, a)
rho = frho(p, T)
Ve = fVe(Vt, rho)
Vetilde = fVetilde(Ve, Ws, W) #W

CN = 2 * W / (rho * Vt**2 * S)

Deltaxcg = fDeltaxcg()
ddeltae_dalpha, uu_intercept, uu_r_value, uu_p_value, uu_std_err = stats.linregress(alpha_rad,CL) # Lots of unused (uu_) values
Cmdelta = - CN * Deltaxcg / Deltadeltae_rad / c
Cmalpha = - ddeltae_dalpha * Cmdelta

Ttotal = fTtotal(T,n_test, hp, M, FFl, FFr)
Ts = fTtotal([T0], 1, [0], [0], mdotfs, mdotfs)
Tcs = Ts * 2 / (rho * Vt**2 * S)
Tc = Ttotal * 2 / (rho * Vt**2 * S)
deltastareeq = delatee_rad - CmTc * (Tcs - Tc) / Cmdelta
Fstareaer = Fe * Ws / W
