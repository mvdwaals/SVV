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

hp = hp_ft * ft_m
Vc = Vc_kts * kts_ms
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

Ttotal = fTtotal(n_test, hp, M, FFl, FFr)
Ts = fTtotal(1, [0], 
