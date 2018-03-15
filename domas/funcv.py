import numpy as np
from dataconst import *

def fp(hp): return p0 * (1 + labda * hp / T0) ** (-g / labda / R)

def fM(p, Vc):
    M_calc_1 = 1 + (gamma - 1) / (2 * gamma) * rho0 / p0 * Vc**2 #()
    M_calc_2 = M_calc_1 ** (gamma / (gamma - 1)) - 1 #{}
    M_calc_3 = 1 + p0 / p * M_calc_2 #()
    M_calc_4 = M_calc_3 ** ((gamma - 1) / gamma) - 1#[]
    M = np.sqrt(2 / (gamma - 1) * M_calc_4)
    return M

def fT(M, Tm): return Tm / (1 + (gamma - 1) / 2 * M**2)

def fa(T): return np.sqrt(gamma * R * T)

def fVt(M, a): return M * a

def frho(p, T): return p / R / T

def fVe(Vt, rho): return Vt * np.sqrt(rho / rho0)

def fVetilde(Ve, Ws, W): return Ve * np.sqrt(Ws / W)
 
