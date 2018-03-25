clear all;
close all;

load('flightdata.mat');
Symmetric_motion_copy_21318;

%Initial time: Dutch Roll - Undamped
dutch_roll_undamped = struct('Name','Dutch_Roll_Undamped','t_init',36550,'t_end',36750,'input',2);

%Initial time: Dutch Roll - Yaw Damper
dutch_roll_damped = struct('Name','Dutch_Roll_Damped','t_init',37185,'t_end',37385,'input',2);

%Initial time: Aperiodic Roll
%aper_roll = struct('Name','Aperiodic Roll','t_init',35555,'t_end',35625);
aper_roll = struct('Name','Aperiodic_Roll','t_init',35569,'t_end',35629,'input',1);

%Initial time: Spiral - 1st
spiral_1 = struct('Name','Spiral_1st','t_init',38250,'t_end',38850,'input',2);

%Initial time: Spiral - 2nd
%spiral_2 = struct('Name','Spiral (2nd)','t_init',38900,'t_end',39500);
spiral_2 = struct('Name','Spiral_2nd','t_init',38800,'t_end',39600,'input',2);

%maneuvers = [dutch_roll_undamped, dutch_roll_damped, spiral_1, spiral_2, aper_roll];
maneuvers = [dutch_roll_undamped, dutch_roll_damped, spiral_2, aper_roll];
%maneuvers = [dutch_roll_undamped];
colorz = get(gca, 'colororder');

for maneuver = maneuvers
    t = maneuver.t_init;
    t_init = maneuver.t_init;
    t_end = maneuver.t_end;

Asymmetric_motion;

V=V0;

C1_asym = [(CYbdot-2*mub)*b/V    	0       0       0;...
            0                   -1/2*b/V	0	0;...
            0                   0           -4*mub*KX2*b^2/(2*V^2)	 4*mub*KXZ*b^2/(2*V^2);...
            Cnbdot*b/V          0            4*mub*KXZ*b^2/(2*V^2)  -4*mub*KZ2*b^2/(2*V^2)];
        
C2_asym = [CYb      CL          CYp*b/(2*V)        (CYr-4*mub)*b/(2*V);...
            0       0           b/(2*V)             0;...
            Clb     0           Clp*b/(2*V)         Clr*b/(2*V);...
            Cnb     0           Cnp*b/(2*V)         Cnr*b/(2*V)];

C3_asym = [CYda  CYdr;...
            0   0;...
            Clda    Cldr;...
            Cnda    Cndr];
        
        
 A = -inv(C1_asym)*C2_asym;
 B = -inv(C1_asym)*C3_asym;
 
 C = eye(4);
 D = zeros(4,2);
 
t = flightdata.time.data(t_init:t_end);
 
 u_input = [(flightdata.delta_a.data(t_init:t_end)'-flightdata.delta_a.data(t_init));...
            -(flightdata.delta_r.data(t_init:t_end)'-flightdata.delta_r.data(t_init))];
        
 display_state = [flightdata.display_active_screen.data(t_init:t_end)';...
                    flightdata.display_graph_state.data(t_init:t_end)'];

 t= t-t(1);
 
 sys = ss(A,B,C,D);
 
testttt = lsim(sys,u_input,t);

filename = maneuver.Name;

%P - R plots
figure();
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 10 14];

%Position definition: left bottom width height
% ax(1) = subplot('Position',[0.125 0.77 0.85 0.22]);
% ax(2) = subplot('Position',[0.125 0.66 0.85 0.10]);
% ax(3) = subplot('Position',[0.125 0.41 0.85 0.22]);
% ax(4) = subplot('Position',[0.125 0.3 0.85 0.10]);
% ax(5) = subplot('Position',[0.125 0.1 0.85 0.16]);

ax(1) = subplot('Position',[0.125 0.77 0.85 0.21]);
flight_data_series = deg2rad(flightdata.Ahrs1_bRollRate.data(t_init:t_end));
model_data_series = deg2rad(testttt(:,3)+flightdata.Ahrs1_bRollRate.data(t_init));
plot(t, flight_data_series, 'Color','r'); hold on;
plot(t, model_data_series,'Color','b'); 
%xlabel('Time [s]');
grid on;
ylim_min = min([min(flight_data_series) min(model_data_series)]);
ylim_max = max([max(flight_data_series) max(model_data_series)]);
ylim_range = ylim_max - ylim_min;
ylim([ylim_min-0.1*ylim_range ylim_max+0.1*ylim_range]);
set(gca,'XTick',[]);
ylabel('p [rad/s]');
hold off;

delta_data_series = model_data_series - flight_data_series;
ax(2) = subplot('Position',[0.125 0.66 0.85 0.10]);
plot(t, delta_data_series,'Color',colorz(5,:)); hold on;
%plot_title = ['Yaw Rate over time - ', num2str(maneuver.Name)];
%title(plot_title);
grid on;
ylim_min = min(delta_data_series);
ylim_max = max(delta_data_series);
ylim_range = ylim_max - ylim_min;
ylim([ylim_min-0.1*ylim_range ylim_max+0.1*ylim_range]);
%xlabel('Time [s]');
ylabel('\Delta{}p [rad/s]');
%legend('Delta - Roll rate','Location','southeast');

model_data_series = deg2rad(testttt(:,4)+flightdata.Ahrs1_bYawRate.data(t_init));
flight_data_series = deg2rad(flightdata.Ahrs1_bYawRate.data(t_init:t_end));
ax(3) = subplot('Position',[0.125 0.41 0.85 0.21]);
mod_plot = plot(t, model_data_series, 'Color','b'); hold on;
flt_plot = plot(t, flight_data_series,'Color','r'); 
ylim_min = min([min(flight_data_series) min(model_data_series)]);
ylim_max = max([max(flight_data_series) max(model_data_series)]);
ylim_range = ylim_max - ylim_min;
ylim([ylim_min-0.1*ylim_range ylim_max+0.1*ylim_range]);
xlabel('Time [s]');
ylabel('r [rad/s]');
set(gca,'XTick',[]);
grid on;
%legend('Simulated Data','Flight Data','Location','southeast'); hold off;

delta_data_series = model_data_series - flight_data_series;
ax(4) = subplot('Position',[0.125 0.3 0.85 0.10]);
del_plot = plot(t, delta_data_series,'Color',colorz(5,:)); hold on;
ylim_min = min(delta_data_series);
ylim_max = max(delta_data_series);
ylim_range = ylim_max - ylim_min;
ylim([ylim_min-0.1*ylim_range ylim_max+0.1*ylim_range]);
%xlabel('Time [s]');
ylabel('\Delta{}r [rad/s]');
%legend('Delta - Yaw rate','Location','southeast');
grid on;

if maneuver.input == 1
    y_ax_inp = ['\delta_a [deg]'];
end
if maneuver.input == 2
    y_ax_inp = ['\delta_r [deg]'];
end

input_data_series = u_input(maneuver.input,:);
ax(5) = subplot('Position',[0.125 0.1 0.85 0.16]);
inp_plot = plot(t, u_input(maneuver.input,:), 'Color','k'); 
% plot_title = ['Yaw Rate over time - ', num2str(maneuver.Name)];
% title(plot_title);
ylim_min = min(input_data_series);
ylim_max = max(input_data_series);
ylim_range = ylim_max - ylim_min;
ylim([ylim_min-0.1*ylim_range ylim_max+0.1*ylim_range]);
% xlabel('Time [s]','Units','normalized','Position',[0.08 -0.25]);
xlabel('Time [s]');
ylabel(y_ax_inp);
%legend('Rudder input','Location','northeast'); 
grid on;

hL = legend([mod_plot,flt_plot, del_plot,inp_plot],{'Simulation Model','Flight Data','Input','Delta'},'Orientation','horizontal');
newPosition = [0.5 0.01 0.01 0.03];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);
saveas(gcf,[filename,'_p_r.png']); 

figure();
fig = gcf;
fig.PaperUnits = 'centimeters';
fig.PaperPosition = [0 0 10 14];
%PHI - R plots


flight_data_series = flightdata.Ahrs1_Roll.data(t_init:t_end);
model_data_series = testttt(:,2)+flightdata.Ahrs1_Roll.data(t_init);
ax(1) = subplot('Position',[0.125 0.77 0.85 0.21]);
plot(t, testttt(:,2)+flightdata.Ahrs1_Roll.data(t_init), 'Color','b'); hold on;
plot(t, flight_data_series,'Color','r'); 
%plot_title = ['Yaw Rate over time - ', num2str(maneuver.Name)];
%title(plot_title);
%xlabel('Time [s]');
grid on;
ylim_min = min([min(flight_data_series) min(model_data_series)]);
ylim_max = max([max(flight_data_series) max(model_data_series)]);
ylim_range = ylim_max - ylim_min;
ylim([ylim_min-0.1*ylim_range ylim_max+0.1*ylim_range]);
ylabel('\phi [deg]');
set(gca,'XTick',[]);
hold off;
%legend('Simulated Data','Flight Data','Location','southeast'); hold off;

delta_data_series = testttt(:,2)+flightdata.Ahrs1_Roll.data(t_init)-flightdata.Ahrs1_Roll.data(t_init:t_end);
ax(2) = subplot('Position',[0.125 0.66 0.85 0.10]);
plot(t, delta_data_series,'Color',colorz(5,:)); hold on;
%plot_title = ['Yaw Rate over time - ', num2str(maneuver.Name)];
%title(plot_title);
ylim_min = min(delta_data_series);
ylim_max = max(delta_data_series);
ylim_range = ylim_max - ylim_min;
ylim([ylim_min-0.1*ylim_range ylim_max+0.1*ylim_range]);
%xlabel('Time [s]');
ylabel('\Delta{}\phi [deg/s]');
%legend('Delta - Bank angle','Location','southeast');
grid on;


model_data_series = deg2rad(testttt(:,4)+flightdata.Ahrs1_bYawRate.data(t_init));
flight_data_series = deg2rad(flightdata.Ahrs1_bYawRate.data(t_init:t_end));
ax(3) = subplot('Position',[0.125 0.41 0.85 0.21]);
mod_plot = plot(t, model_data_series, 'Color','b'); hold on;
flt_plot = plot(t, flight_data_series,'Color','r'); 
ylim_min = min([min(flight_data_series) min(model_data_series)]);
ylim_max = max([max(flight_data_series) max(model_data_series)]);
ylim_range = ylim_max - ylim_min;
ylim([ylim_min-0.1*ylim_range ylim_max+0.1*ylim_range]);
grid on;
%xlabel('Time [s]');
ylabel('r [rad/s]');
set(gca,'XTick',[]);
%legend('Simulated Data','Flight Data','Location','southeast'); hold off;

delta_data_series = deg2rad(testttt(:,4)+flightdata.Ahrs1_bYawRate.data(t_init)-flightdata.Ahrs1_bYawRate.data(t_init:t_end));
ax(4) = subplot('Position',[0.125 0.3 0.85 0.10]);
del_plot = plot(t, delta_data_series,'Color',colorz(5,:)); hold on;
ylim_min = min(delta_data_series);
ylim_max = max(delta_data_series);
ylim_range = ylim_max - ylim_min;
ylim([ylim_min-0.1*ylim_range ylim_max+0.1*ylim_range]);
%xlabel('Time [s]','Units','normalized','Position',[0.07 -0.2]);
ylabel('\Delta{}r [rad/s]');
%legend('Delta - Yaw rate','Location','southeast');
grid on;

if maneuver.input == 1
    y_ax_inp = ['\delta_a [deg]'];
end
if maneuver.input == 2
    y_ax_inp = ['\delta_r [deg]'];
end

input_data_series = u_input(maneuver.input,:);
ax(5) = subplot('Position',[0.125 0.1 0.85 0.16]);
inp_plot = plot(t, u_input(maneuver.input,:), 'Color','k'); 
% plot_title = ['Yaw Rate over time - ', num2str(maneuver.Name)];
% title(plot_title);
ylim_min = min(input_data_series);
ylim_max = max(input_data_series);
ylim_range = ylim_max - ylim_min;
ylim([ylim_min-0.1*ylim_range ylim_max+0.1*ylim_range]);
%xlabel('Time [s]','Units','normalized','Position',[0.08 -0.25]);
xlabel('Time [s]');
ylabel(y_ax_inp);
%legend('Rudder input','Location','northeast'); 
grid on; 

hL = legend([mod_plot,flt_plot, del_plot,inp_plot],{'Simulation Model','Flight Data','Input','Delta'},'Orientation','horizontal');
newPosition = [0.5 0.01 0.01 0.03];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits);

saveas(gcf,[filename,'_phi_r.png']);
end

opt = pzoptions;
opt.Ylim = [-3 3];

figure('Visible','on');
pzplot(sys);