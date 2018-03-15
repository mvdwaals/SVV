import numpy as np

# Added
p0 = 101325 #[N/m2]
pi = np.pi
gamma = 1.4
Ws = 90500 #N

ft_m = 0.3048
kts_ms = 0.514444
c_k = 273.15
lbhr_kgs = 0.000125998
lb_kg = 0.453592
in_m = 0.0254
mdotfs = 0.048 #kg/s
CmTc = -0.0064

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
T0  = 288.15          # temperature at sea level in ISA [K]
R      = 287.05          # specific gas constant [m^2/sec^2K]
g      = 9.81            # [m/sec^2] (gravity constant)
