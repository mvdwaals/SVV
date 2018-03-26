%% CALCULATING AIRCRAFT MASS AS A FUNCTION OF TIME
% Calculating the initial mass
m_extra_payload = convmass(220, 'lbm', 'kg');
m_fuel_start = convmass(2650, 'lbm', 'kg');
m_people = 92 + 95 + 76 + 61 + 59 + 66 + 77 + 77 + 84;
m_OEW = convmass(9165, 'lbm', 'kg');

m_start = m_OEW + m_people + m_extra_payload + m_fuel_start;

% Calculate the decrease in mass as a function of time
fuelusedl_kg = convmass (flightdata.lh_engine_FU.data,'lbm', 'kg');
fuelusedr_kg = convmass (flightdata.rh_engine_FU.data, 'lbm', 'kg');

ac_mass_lst = m_start - fuelusedl_kg - fuelusedr_kg;