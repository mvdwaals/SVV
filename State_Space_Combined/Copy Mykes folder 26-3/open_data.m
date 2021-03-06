alpha = flightdata.vane_AOA.data;
t = flightdata.time.data;
dte = flightdata.elevator_dte.data;
Fe = flightdata.column_fe.data;
FFl = flightdata.lh_engine_FMF.data;
FFr = flightdata.rh_engine_FMF.data;
fuelusedl = flightdata.lh_engine_FU.data;
fuelusedr = flightdata.rh_engine_FU.data;
delta_a = flightdata.delta_a.data;
delta_e = flightdata.delta_e.data;
delta_r = flightdata.delta_r.data;
roll = flightdata.Ahrs1_Roll.data;
pitch = flightdata.Ahrs1_Pitch.data;
rollrate = flightdata.Ahrs1_bRollRate.data;
pitchrate = flightdata.Ahrs1_bPitchRate.data;
yawrate = flightdata.Ahrs1_bYawRate.data;
longacc = flightdata.Ahrs1_bLongAcc.data;
latacc = flightdata.Ahrs1_bLatAcc.data;
normacc = flightdata.Ahrs1_bNormAcc.data;
vertacc = flightdata.Ahrs1_VertAcc.data;
ahdgacc = flightdata.Ahrs1_aHdgAcc.data;
xhdgacc = flightdata.Ahrs1_xHdgAcc.data;
sat = flightdata.Dadc1_sat.data;
tat = flightdata.Dadc1_tat.data;
h = flightdata.Dadc1_alt.data;
hbc = flightdata.Dadc1_bcAlt.data;
mach = flightdata.Dadc1_mach.data;
Vcas = flightdata.Dadc1_cas.data;
Vtas = flightdata.Dadc1_tas.data;
hrate = flightdata.Dadc1_altRate.data;
measure = flightdata.measurement_running.data;

%measuremtns: 1060 - 1900
%short period: 3200 - 3250
%phugoid: 3250 - 3500
%aperiodic roll: 3550 - 3600
%Dutch roll without damper: 3660 - 3690
%Dutch roll with damper: 3725 - 3745
%spiral: 3890 - 4050

tx = t(10511:18911); %unit?
hx = h(10511:18911); % ft
Vcasx = Vcas(10511:18911); % kts
tatx = tat(10511:18911); % deg Celsius
FFlx = FFl(10511:18911); % lbs/hr
FFrx = FFr(10511:18911); % lbs/hr
fuelusedlx = fuelusedl(10511:18911); % lbs
fuelusedrx = fuelusedr(10511:18911); % lbs
alphax = alpha(10511:18911); % deg







