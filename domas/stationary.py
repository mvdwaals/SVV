import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import os 

pi = np.pi
gamma = 1.4

ft_m = 0.3048
kts_ms = 0.514444
c_k = 273.15
lbshr_kgs = 0.000125998
lbs_kg = 0.453592

# h[ft] v[kts] tat[c] ffl[lbs/hr] ffr[lbs/hr] wf[lbs] alpha[deg]
data_not_si = np.array([#[5010, 248, 11.2, 768, 784, 215, 1.5],
                        [6040, 251, 10.5, 786, 776, 410, 1.3],
                        [6035, 222, 8.2, 637, 650, 450, 2.0],
                        [6030, 192, 6.1, 515, 548, 487, 3.0],
                        [6040, 165, 4.8, 558, 470, 540, 4.5],
                        [6030, 130, 3.5, 373, 385, 574, 8.3],
                        [6040, 115, 3.0, 414, 435, 604, 10.1]])

n_tests = len(data_not_si)

data_not_si_T = data_not_si.T
h_ft = data_not_si_T[0]
V_kts = data_not_si_T[1]
TAT_c = data_not_si_T[2]
FFl_lbshr = data_not_si_T[3]
FFr_lbshr = data_not_si_T[4]
Wf_lbs = data_not_si_T[5]
alpha_deg = data_not_si_T[6]

# h[m] v[m/s] tat[K] ffl[kg/s] ffr[kg/s] wf[kg] alpha[rad]
person_weight_kg = [92, 95, 76, 61, 59, 66, 77, 77, 84]
empty_weight_lbs = 9165
total_fuel_lbs = 2650

empty_weight_kg = empty_weight_lbs * lbs_kg
total_fuel_kg = total_fuel_lbs * lbs_kg

m_tot = sum(person_weight_kg) + empty_weight_kg + total_fuel_kg

hp0 = h_ft * ft_m
V0 = V_kts * kts_ms
alpha0 = np.radians(alpha_deg)
th0 = alpha0 # Slope was zero for all these measurements
m = m_tot - Wf_lbs * lbs_kg
Tmeas = TAT_c + c_k
FFl = FFl_lbshr * lbshr_kgs
FFr = FFr_lbshr * lbshr_kgs

# xcg = 0.25 * c    # Not sure what to do with this

# Aircraft geometry
S      = 30.00	          # wing area [m^2]
Sh     = 0.2 * S         # stabiliser area [m^2]
Sh_S   = Sh / S	          # [ ]
lh     = 0.71 * 5.968    # tail length [m]
c      = 2.0569	          # mean aerodynamic cord [m]
lh_c   = lh / c	          # [ ]
b      = 15.911	          # wing span [m]
bh     = 5.791	          # stabilser span [m]
A      = b ** 2 / S      # wing aspect ratio [ ]
Ah     = bh ** 2 / Sh    # stabilser aspect ratio [ ]
Vh_V   = 1	          # [ ]
ih     = -2 * pi / 180   # stabiliser angle of incidence [rad]

# Constant values concerning atmosphere and gravity
rho0   = 1.2250          # air density at sea level [kg/m^3] 
labda = -0.0065         # temperature gradient in ISA [K/m]
Temp0  = 288.15          # temperature at sea level in ISA [K]
R      = 287.05          # specific gas constant [m^2/sec^2K]
g      = 9.81            # [m/sec^2] (gravity constant)

# Added
p0 = 101325 #[N/m2]
a = 343 #[m/s]

# air density [kg/m^3]  
rho    = rho0 * ((1+(labda * hp0 / Temp0)))**(-((g / (labda*R)) + 1))   
W      = m * g            # [N]       (aircraft weight)

# ISA calculations
p = p0 * (1 + labda * hp0 / Temp0) ** (-g / labda / R)

M_calc_1 = 1 + (gamma - 1) / (2 * gamma) * rho0 / p0 * V0**2
M_calc_2 = M_calc_1**(gamma / (gamma - 1)) - 1
M_calc_3 = (1 + p0/p * M_calc_2) ** ((gamma - 1) / gamma)
M = np.sqrt(2 / (gamma - 1) * (M_calc_3 - 1))

T = Tmeas / (1 + (gamma - 1) / 2 * M**2)

Vt = M * a
Ve = Vt * np.sqrt(rho / rho0)

# Constant values concerning aircraft inertia

muc    = m / (rho * S * c)
mub    = m / (rho * S * b)
KX2    = 0.019
KZ2    = 0.042
KXZ    = 0.002
KY2    = 1.25 * 1.114

# Lift coefficient
CL = 2 * W / (rho * V0 ** 2 * S)              # Lift coefficient [ ]

#Plotting and finding CLa
CLa, intercept, uu_r_value, uu_p_value, uu_std_err = stats.linregress(alpha_deg,CL) # Lots of unused (uu_) values

linregress_x = np.array([min(alpha_deg), max(alpha_deg)])
linregress_y = intercept + CLa * linregress_x

plt.plot(linregress_x, linregress_y)
plt.scatter(alpha_deg, CL)
plt.ylabel('CL [-]')
plt.xlabel('alpha [deg]')
plt.savefig('CL_alphadeg.png')

# CD calculations
T_ISA = Temp0 + labda * hp0
DeltaT = abs(T_ISA - T)

outfile = open('matlab.dat', 'w')
for i in range(n_tests):
    outfile.write(str(hp0[i])+' '+str(M[i])+' '+str(DeltaT[i])+' '+str(FFl[i])+' '+str(FFr[i])+'\n')
outfile.close()

# This doesn't work on my device, sadly
# os.system('java -jar thrust.jar')

dummy = input("Only continue if thrust.dat exists.")

infile = np.genfromtxt('thrust.dat').T
T1 = infile[0]
T2 = infile[1]
Ttotal = T1 + T2

CD = CD0 + (CLa * alpha0) ** 2 / (pi * A * e) # Drag coefficient [ ]

# Aerodynamic constants

Cmac   = 0                      # Moment coefficient about the aerodynamic centre [ ]
CNwa   = CLa                    # Wing normal force slope [ ]
CNha   = 2 * pi * Ah / (Ah + 2) # Stabiliser normal force slope [ ]
depsda = 4 / (A + 2)            # Downwash gradient [ ]

# Stabiblity derivatives

CX0    = W * sin(th0) / (0.5 * rho * V0 ** 2 * S)
CXu    = -0.02792
CXa    = -0.47966
CXadot = +0.08330
CXq    = -0.28170
CXde   = -0.03728

CZ0    = -W * cos(th0) / (0.5 * rho * V0 ** 2 * S)
CZu    = -0.37616
CZa    = -5.74340
CZadot = -0.00350
CZq    = -5.66290
CZde   = -0.69612

Cmu    = +0.06990
Cmadot = +0.17800
Cmq    = -8.79415

CYb    = -0.7500
CYbdot =  0     
CYp    = -0.0304
CYr    = +0.8495
CYda   = -0.0400
CYdr   = +0.2300

Clb    = -0.10260
Clp    = -0.71085
Clr    = +0.23760
Clda   = -0.23088
Cldr   = +0.03440

Cnb    =  +0.1348
Cnbdot =   0     
Cnp    =  -0.0602
Cnr    =  -0.2061
Cnda   =  -0.0120
Cndr   =  -0.0939

