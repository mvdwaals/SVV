clear all;
close all;

load('flightdata.mat');
Symmetric_motion_copy_21318;

%Initial time: Dutch Roll - Undamped
dutch_roll_undamped = struct('Name','Dutch Roll - Undamped','t_init',36550,'t_end',36750,'input',2);

%Initial time: Dutch Roll - Yaw Damper
dutch_roll_damped = struct('Name','Dutch Roll - Damped','t_init',37185,'t_end',37335,'input',2);

%Initial time: Aperiodic Roll
%aper_roll = struct('Name','Aperiodic Roll','t_init',35555,'t_end',35625);
aper_roll = struct('Name','Aperiodic Roll','t_init',35569,'t_end',35625,'input',1);

%Initial time: Spiral - 1st
spiral_1 = struct('Name','Spiral (1st)','t_init',38250,'t_end',38850,'input',2);

%Initial time: Spiral - 2nd
%spiral_2 = struct('Name','Spiral (2nd)','t_init',38900,'t_end',39500);
spiral_2 = struct('Name','Spiral (2nd)','t_init',38800,'t_end',39600,'input',2);

%maneuvers = [dutch_roll_undamped, dutch_roll_damped, spiral_1, spiral_2, aper_roll];
maneuvers = [dutch_roll_undamped, dutch_roll_damped, spiral_2, aper_roll];
%maneuvers = [dutch_roll_undamped];

for maneuver = maneuvers
    t = maneuver.t_init;
    t_init = maneuver.t_init;
    t_end = maneuver.t_end;

%maneuver = struct(
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

figure();
ax(1) = subplot(3,1,1);
plot(t, deg2rad(testttt(:,3)+flightdata.Ahrs1_bRollRate.data(t_init)), 'Color','r'); hold on;
plot(t, deg2rad(flightdata.Ahrs1_bRollRate.data(t_init:t_end)),'Color','b'); 
%plot_title = ['Yaw Rate over time - ', num2str(maneuver.Name)];
%title(plot_title);
xlabel('Time [s]');
ylabel('p [rad/s]');
grid on;
legend('Simulated Data','Flight Data','Location','southeast'); hold off;

ax(2) = subplot(3,1,2);
plot(t, abs(testttt(:,3)+flightdata.Ahrs1_bRollRate.data(t_init)-deg2rad(flightdata.Ahrs1_bRollRate.data(t_init:t_end))),'Color','b'); hold on;
%plot_title = ['Yaw Rate over time - ', num2str(maneuver.Name)];
%title(plot_title);
xlabel('Time [s]');
ylabel('p [rad/s]');
legend('Delta - Roll rate','Location','southeast');
grid on;

ax(3) = subplot(3,1,3);
plot(t, u_input(maneuver.input,:), 'Color','k'); hold on;
% plot_title = ['Yaw Rate over time - ', num2str(maneuver.Name)];
% title(plot_title);
xlabel('Time [s]');
ylabel('\delta_a [deg]');
legend('Aileron input','Location','southeast'); 
grid on; 

linkaxes(ax,'x')
hold off;

saveas(gcf,[filename,'_ROLL.png']);

figure();
ax(1) = subplot(3,1,1);
plot(t, deg2rad(testttt(:,4)+flightdata.Ahrs1_bYawRate.data(t_init)), 'Color','r'); hold on;
plot(t, deg2rad(flightdata.Ahrs1_bYawRate.data(t_init:t_end)),'Color','b'); 
%plot_title = ['Yaw Rate over time - ', num2str(maneuver.Name)];
%title(plot_title);
xlabel('Time [s]');
ylabel('r [rad/s]');
grid on;
legend('Simulated Data','Flight Data','Location','southeast'); hold off;

ax(2) = subplot(3,1,2);
plot(t, abs(deg2rad(testttt(:,4)+flightdata.Ahrs1_bYawRate.data(t_init)-flightdata.Ahrs1_bYawRate.data(t_init:t_end))),'Color','b'); hold on;
%plot_title = ['Yaw Rate over time - ', num2str(maneuver.Name)];
%title(plot_title);
xlabel('Time [s]');
ylabel('r [rad/s]');
legend('Delta - Yaw rate','Location','southeast');
grid on;

ax(3) = subplot(3,1,3);
plot(t, u_input(maneuver.input,:), 'Color','k'); hold on;
% plot_title = ['Yaw Rate over time - ', num2str(maneuver.Name)];
% title(plot_title);
xlabel('Time [s]');
ylabel('\delta_a [deg]');
legend('Aileron input','Location','southeast'); 
grid on; 

linkaxes(ax,'x');

saveas(gcf,[filename,'_YAW.png']);

hold off;
end

opt = pzoptions;
opt.Ylim = [-3 3];

figure('Visible','on');
pzplot(sys);