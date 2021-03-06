clear all;
close all;

%Load flightdata
load('flightdata.mat');

%Initialize maneuvers
maneuver_initialize;

%Load Matlab Default colors for plots
palette = get(gca, 'colororder');

%Load aircraft mass as a function of time
ac_mass_time;

for maneuver = maneuvers
t_init = maneuver.t_init;
t_end = maneuver.t_end;

%Load stability derivatives and data from stationary measurements
Cit_par;

V=V0;

if strcmp(maneuver.mode,['a'])
    Matrices_a;
     u_input = [(flightdata.delta_a.data(t_init:t_end)'-flightdata.delta_a.data(t_init));...
            -(flightdata.delta_r.data(t_init:t_end)'-flightdata.delta_r.data(t_init))];
end

if strcmp(maneuver.mode,['s'])
    Matrices_s;
     u_input = [deg2rad(flightdata.delta_e.data(t_init:t_end))'-deg2rad(flightdata.delta_e.data(t_init))];
end

%Load maneuver time axes and set initial time to zero
t = flightdata.time.data(t_init:t_end);
t= t-t(1);
 
%Initialize state space system
sys = ss(A,B,C,D);
opt = pzoptions;
opt.Ylim = [-3 3];

figure(3);
pzplot(sys);

%Calculate response
response = lsim(sys,u_input,t);

filename = maneuver.Name;

%Plotting commences
figure();
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 10 14];

if strcmp(maneuver.first_plot,'p_roll_rate')
    flight_data_series1 = flightdata.Ahrs1_bRollRate.data(t_init:t_end);
    model_data_series1 = response(:,3)+flightdata.Ahrs1_bRollRate.data(t_init);
    delta_data_series1 = model_data_series1 - flight_data_series1;
    y_label_1 = 'p [deg/s]';
    y_label_1_delta = '\Delta{}p [deg/s]';
end

if strcmp(maneuver.first_plot,'phi_bank_angle')
    flight_data_series1 = flightdata.Ahrs1_Roll.data(t_init:t_end);
    model_data_series1 = response(:,2)+flightdata.Ahrs1_Roll.data(t_init);
    delta_data_series1 = model_data_series1 - flight_data_series1;
    y_label_1 = '\phi [deg]';
    y_label_1_delta = '\Delta{}\phi [deg]';
end

if strcmp(maneuver.first_plot,'theta_pitch_angle')
    flight_data_series1 = flightdata.Ahrs1_Pitch.data(t_init:t_end);
    model_data_series1 = rad2deg(response(:,3))+flightdata.Ahrs1_Pitch.data(t_init);
    delta_data_series1 = model_data_series1 - flight_data_series1;
    y_label_1 = '\theta [deg]';
    y_label_1_delta = '\Delta{}\theta [deg]';
end

if strcmp(maneuver.second_plot,'r_yaw_rate')
    model_data_series2 = response(:,4)+flightdata.Ahrs1_bYawRate.data(t_init);
    flight_data_series2 = flightdata.Ahrs1_bYawRate.data(t_init:t_end);
    delta_data_series2 = model_data_series2 - flight_data_series2;
    y_label_2 = 'r [deg/s]';
    y_label_2_delta = '\Delta{}r [deg/s]';
end

if strcmp(maneuver.second_plot,'V')
    model_data_series2 = response(:,1)+flightdata.Dadc1_tas.data(t_init)*0.51444;
    flight_data_series2 = flightdata.Dadc1_tas.data(t_init:t_end)*0.51444;
    delta_data_series2 = model_data_series2 - flight_data_series2;
    y_label_2 = 'V [m/s]';
    y_label_2_delta = '\Delta{}V [m/s]';
end

if strcmp(maneuver.second_plot,'q_pitch_rate')
    model_data_series2 = rad2deg(response(:,4))+flightdata.Ahrs1_bPitchRate.data(t_init);
    flight_data_series2 = flightdata.Ahrs1_bPitchRate.data(t_init:t_end);
    delta_data_series2 = model_data_series2 - flight_data_series2;
    y_label_2 = 'q [deg/s]';
    y_label_2_delta = '\Delta{}q [deg/s]';
end

if strcmp(maneuver.second_plot,'p_roll_rate')
    flight_data_series2 = flightdata.Ahrs1_bRollRate.data(t_init:t_end);
    model_data_series2 = response(:,3)+flightdata.Ahrs1_bRollRate.data(t_init);
    delta_data_series2 = model_data_series1 - flight_data_series1;
    y_label_2 = 'p [deg/s]';
    y_label_2_delta = '\Delta{}p [deg/s]';
end

%Load input


if maneuver.input == 1
    y_ax_inp = ['\delta_a [deg]'];
    input_data_series = u_input(maneuver.input,:);
end
if maneuver.input == 2
    y_ax_inp = ['\delta_r [deg]'];
    input_data_series = u_input(maneuver.input,:);
end
if maneuver.input == 3
    y_ax_inp = ['\delta_e [deg]'];
    maneuver.input = 1;
    input_data_series = u_input(maneuver.input,:);
end



ax(1) = subplot('Position',[0.125 0.77 0.85 0.21]);
plot(t, flight_data_series1, 'Color','r'); hold on;
plot(t, model_data_series1,'Color','b'); 
grid on;
ylim_min = min([min(flight_data_series1) min(model_data_series1)]);
ylim_max = max([max(flight_data_series1) max(model_data_series1)]);
ylim_range = ylim_max - ylim_min;
ylim([ylim_min-0.1*ylim_range ylim_max+0.1*ylim_range]);
set(gca,'xticklabel',[]);
ylabel(y_label_1);
hold off;

ax(2) = subplot('Position',[0.125 0.66 0.85 0.10]);
plot(t, delta_data_series1,'Color',palette(5,:)); hold on;
grid on;
ylim_min = min(delta_data_series1);
ylim_max = max(delta_data_series1);
ylim_range = ylim_max - ylim_min;
ylim([ylim_min-0.1*ylim_range ylim_max+0.1*ylim_range]);
ylabel(y_label_1_delta);

ax(3) = subplot('Position',[0.125 0.41 0.85 0.21]);
mod_plot = plot(t, model_data_series2, 'Color','b'); hold on;
flt_plot = plot(t, flight_data_series2,'Color','r'); 
ylim_min = min([min(flight_data_series2) min(model_data_series2)]);
ylim_max = max([max(flight_data_series2) max(model_data_series2)]);
ylim_range = ylim_max - ylim_min;
ylim([ylim_min-0.1*ylim_range ylim_max+0.1*ylim_range]);
xlabel('Time [s]');
ylabel(y_label_2);
set(gca,'xticklabel',[]);
grid on;

ax(4) = subplot('Position',[0.125 0.3 0.85 0.10]);
del_plot = plot(t, delta_data_series2,'Color',palette(5,:)); hold on;
ylim_min = min(delta_data_series2);
ylim_max = max(delta_data_series2);
ylim_range = ylim_max - ylim_min;
ylim([ylim_min-0.1*ylim_range ylim_max+0.1*ylim_range]);
ylabel(y_label_2_delta);
grid on;

ax(5) = subplot('Position',[0.125 0.1 0.85 0.16]);
inp_plot = plot(t, input_data_series, 'Color','k'); 
ylim_min = min(input_data_series);
ylim_max = max(input_data_series);
ylim_range = ylim_max - ylim_min;
ylim([ylim_min-0.1*ylim_range ylim_max+0.1*ylim_range]);
xlabel('Time [s]');
ylabel(y_ax_inp);
grid on;

hL = legend([mod_plot,flt_plot, del_plot,inp_plot],{'Simulation Model','Flight Data','Delta','Input'},'Orientation','horizontal');
newPosition = [0.5 0.01 0.01 0.03];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
saveas(gcf,[filename,'_',maneuver.first_plot,'_',maneuver.second_plot,'.png']); 
hold off;
end

