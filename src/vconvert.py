from isadata import *
import numpy as np

def hp_to_p(hp): return p0 * (1 + labda * hp / T0) ** -g0 / labda / R

def p_Vc_to_M(p, Vc):
    M_calc_1 = 1 + (gamma - 1) / (2 * gamma) * rho0 / p0 * Vc**2
    M_calc_2 = M_calc_1 ** (gamma / (gamma - 1)) - 1
    M_calc_3 = (1 + p0/p * M_calc_2) ** ((gamma - 1) / gamma) - 1
    M = np.sqrt(2 / (gamma - 1) * M_calc_3)
    return M

def M_Tm_to_T(M, Tm): return Tm / (1 + (gamma - 1) / 2 * M**2)

def T_to_a(T): return np.sqrt(gamma * R * T)

def M_a_to_Vt(M, a): return M * a

def p_T_to_rho(p, T): return p / R / T

def Vt_rho_to_Ve(Vt, rho): return np.sqrt(rho / rho0)
    
def hp_Vc_Tm_to_Ve(hp, Vc, Tm):
    p = hp_to_p(hp)
    M = p_Vc_to_M(p, Vc)
    T = M_Tm_to_T(M, Tm)
    a = T_to_a(T)
    Vt = M_a_to_Vt(M, a)
    rho = p_T_to_rho(p, T)
    return Vt_rho_to_Ve(Vt, rho)
