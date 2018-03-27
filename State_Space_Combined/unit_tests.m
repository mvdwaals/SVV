clc
clearvars
close all

load('flightdata.mat')

%% maneuver_initialize
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
fprintf('|                      maneuver_initialize                    |\n')
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
maneuver_initialize

fprintf('\nVerifying display overlap\n')
fprintf('========================\n')
for maneuver = maneuvers
    i_disp = flightdata.display_active_screen.data == maneuver.disp;
    
    check = zeros(size(i_disp));
    check(maneuver.t_init:maneuver.t_end) = 1;
    
    res = check .* i_disp;
    if sum(res) > 0
        fprintf('[pass]\tCheck display var: %s\n', maneuver.Name)
    else
        fprintf('[fail]\tCheck display var: %s\n', maneuver.Name)
    end
    
    clear i_disp check res
end

fprintf('\nVerifying timing overlap\n')
fprintf('========================\n')

check = zeros(size(flightdata.time.data));

for maneuver = maneuvers
    check(maneuver.t_init:maneuver.t_end) = check(maneuver.t_init:maneuver.t_end) + 1;
end

if find(check>1)
    fprintf('[fail]\tCheck overlap\n')
else
    fprintf('[pass]\tCheck overlap\n')
end
clear check

fprintf('\nVerifying input signal\n')
fprintf('=======================\n')

for maneuver = maneuvers
    da = flightdata.delta_a.data(maneuver.t_init:maneuver.t_end);
    dr = flightdata.delta_r.data(maneuver.t_init:maneuver.t_end);
    de = flightdata.delta_e.data(maneuver.t_init:maneuver.t_end);
    inputs = [ da' - da(1);
               dr' - dr(1);
               de' - de(1)];
    abs_inputs = abs(inputs);
    avg_inputs = mean(abs_inputs,2);
    [~,idx] = max(avg_inputs);
    
    if idx == maneuver.input
        fprintf('[pass]\tCheck input: %s\n', maneuver.Name)
    else
        fprintf('[fail]\tCheck input: %s\n', maneuver.Name)
    end
    
    clear da dr de inputs abs_inputs avg_inputs val idx
end

clear aper_roll dutch_roll_damped dutch_roll_undamped maneuver maneuvers
clear maneuvers_asym maneuvers_sym phugoid short_period spiral_1 spiral_2

%% ac_mass_time
fprintf('\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
fprintf('|                         ac_mass_time                        |\n')
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')

ac_mass_time

fprintf('\nVerifying decreasing mass\n')
fprintf('=========================\n')
if find(diff(ac_mass_lst)>0)
    fprintf('[fail]\tCheck mass decreasing\n')
else
    fprintf('[pass]\tCheck mass decreasing\n')
end

clear ac_mass_lst fuelusedl_kg fuelusedr_kg m_extra_payload m_fuel_start
clear m_OEW m_people m_start

%% matrices
fprintf('\n+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
fprintf('|                             matrices                        |\n')
fprintf('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')

t_init = 36000;
ac_mass_time

Cit_par
V=V0;

fprintf('\nVerifying symetric matrices\n')
fprintf('============================\n')
Matrices_s

[ma,na] = size(A);
[mb,nb] = size(B);
[mc,nc] = size(C);
[md,nd] = size(D);

if ma == 4 && na == 4
    fprintf('[pass]\tShape of A\n')
else
    fprintf('[fail]\tShape of A\n')
end

if mb == 4 && nb == 1
    fprintf('[pass]\tShape of B\n')
else
    fprintf('[fail]\tShape of B\n')
end

if mc == 4 && nc == 4
    fprintf('[pass]\tShape of C\n')
else
    fprintf('[fail]\tShape of C\n')
end

if md == 4 && nd == 1
    fprintf('[pass]\tShape of D\n')
else
    fprintf('[fail]\tShape of D\n')
end

fprintf('\nVerifying asymetric matrices\n')
fprintf('=============================\n')
Matrices_a

[ma,na] = size(A);
[mb,nb] = size(B);
[mc,nc] = size(C);
[md,nd] = size(D);

if ma == 4 && na == 4
    fprintf('[pass]\tShape of A\n')
else
    fprintf('[fail]\tShape of A\n')
end

if mb == 4 && nb == 2
    fprintf('[pass]\tShape of B\n')
else
    fprintf('[fail]\tShape of B\n')
end

if mc == 4 && nc == 4
    fprintf('[pass]\tShape of C\n')
else
    fprintf('[fail]\tShape of C\n')
end

if md == 4 && nd == 2
    fprintf('[pass]\tShape of D\n')
else
    fprintf('[fail]\tShape of D\n')
end
