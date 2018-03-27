clear all;
close all;
clc;

%Load Matlab Default colors for plots
palette = get(gca, 'colororder');

load('flightdata.mat')
ac_mass_time;

t_init = 38800;
t_end = 39600;
t_maneuver = t_init:t_end;

t = flightdata.time.data(t_maneuver) - flightdata.time.data(t_init);

Cit_par;
V = V0;

da = flightdata.delta_a.data(t_maneuver);
dr = flightdata.delta_r.data(t_maneuver);

U = [ (da'-da(1));...
     -(dr'-dr(1))];
       

Matrices_a;
sys = ss(A,B,C,D);
Y_orig = lsim(sys, U, t);


Clb = -0.1326;
Matrices_a;
sys = ss(A,B,C,D);
Y_tuned = lsim(sys, U, t);

roll_meas  = flightdata.Ahrs1_Roll.data(t_maneuver);
roll_orig  = Y_orig(:,2)  + roll_meas(1);
roll_tuned = Y_tuned(:,2) + roll_meas(1);

roll_orig_delta = roll_orig - roll_meas;
roll_tuned_delta = roll_tuned - roll_meas;

r_meas  = flightdata.Ahrs1_bYawRate.data(t_maneuver);
r_orig  = Y_orig(:,4)  + r_meas(1);
r_tuned = Y_tuned(:,4) + r_meas(1);

r_orig_delta  = r_orig - r_meas;
r_tuned_delta = r_tuned - r_meas;

figure();
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 10 14];

ax(1) = subplot('Position',[0.125 0.77 0.85 0.21]);
plot(t, roll_meas,'r'); hold on;
plot(t, roll_orig,'b');
plot(t, roll_tuned,'b--');
grid on;
ylim_min = min([min(roll_meas) min(roll_orig) min(roll_tuned)]);
ylim_max = max([max(roll_meas) max(roll_orig) max(roll_tuned)]);
ylim_range = ylim_max - ylim_min;
ylim([ylim_min-0.1*ylim_range ylim_max+0.1*ylim_range]);
set(gca,'xticklabel',{[]})
ylabel('\phi [deg]');
hold off;

ax(2) = subplot('Position',[0.125 0.66 0.85 0.10]);
plot(t, roll_orig_delta,'Color',palette(5,:)); hold on;
plot(t, roll_tuned_delta,'--','Color',palette(5,:));
grid on;
ylim_min = min([min(roll_orig_delta) min(roll_tuned_delta)]);
ylim_max = max([max(roll_orig_delta) max(roll_tuned_delta)]);
ylim_range = ylim_max - ylim_min;
ylim([ylim_min-0.1*ylim_range ylim_max+0.1*ylim_range]);
ylabel('\Delta{}\phi [deg]');

ax(3) = subplot('Position',[0.125 0.41 0.85 0.21]);
mod_plot = plot(t, r_meas, 'r'); hold on;
flt_plot = plot(t, r_orig, 'b'); 
tun_plot = plot(t, r_tuned, 'b--');
ylim_min = min([min(r_meas) min(r_orig) min(r_tuned)]);
ylim_max = max([max(r_meas) max(r_orig) max(r_tuned)]);
ylim_range = ylim_max - ylim_min;
ylim([ylim_min-0.1*ylim_range ylim_max+0.1*ylim_range]);
xlabel('Time [s]');
ylabel('r [deg/s]');
set(gca,'xticklabel',{[]}) 
grid on;

ax(4) = subplot('Position',[0.125 0.3 0.85 0.10]);
del_plot = plot(t, r_orig_delta,'Color',palette(5,:)); hold on;
plot(t, r_tuned_delta,'--', 'Color',palette(5,:));
ylim_min = min([min(r_orig_delta) min(r_tuned_delta)]);
ylim_max = max([max(r_orig_delta) max(r_tuned_delta)]);
ylim_range = ylim_max - ylim_min;
ylim([ylim_min-0.1*ylim_range ylim_max+0.1*ylim_range]);
ylabel('\Delta{}r [deg/s]');
grid on;

ax(5) = subplot('Position',[0.125 0.1 0.85 0.16]);
inp_plot = plot(t, dr, 'Color','k'); 
ylim_min = min(dr);
ylim_max = max(dr);
ylim_range = ylim_max - ylim_min;
ylim([ylim_min-0.1*ylim_range ylim_max+0.1*ylim_range]);
xlabel('Time [s]');
ylabel('\delta_r [deg]');
grid on;

hL = legend([mod_plot,flt_plot, del_plot,inp_plot],{'Simulation Model','Flight Data','Delta','Input'},'Orientation','horizontal');
newPosition = [0.5 0.01 0.01 0.03];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
saveas(gcf,['spiral_tuing.png']); 
hold off;
