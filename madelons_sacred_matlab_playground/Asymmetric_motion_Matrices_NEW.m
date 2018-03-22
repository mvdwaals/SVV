clear all;

load('flightdata.mat');

%Initial time: Aperiodic Roll
aper_roll = struct('Name','Aperiodic Roll','t_init',35500,'t_end',36100);

%Initial time: Dutch Roll - Undamped
dutch_roll_undamped = struct('Name','Dutch Roll - Undamped','t_init',36520,'t_end',36820);

%Initial time: Dutch Roll - Yaw Damper
dutch_roll_damped = struct('Name','Dutch Roll - Damped','t_init',37150,'t_end',37450);

%Initial time: Spiral - 1st
spiral_1 = struct('Name','Spiral (1st)','t_init',38250,'t_end',38850);

%Initial time: Spiral - 2nd
%spiral_2 = struct('Name','Spiral (2nd)','t_init',38900,'t_end',39500);
spiral_2 = struct('Name','Spiral (2nd)','t_init',38800,'t_end',39500);

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

 t= t-t(1);
 
 sys = ss(A,B,C,D);
 
%  x0 = [0;
%      deg2rad(flightdata.Ahrs1_Roll.data(t_init));
%      flightdata.Ahrs1_bRollRate.data(t_init);
%      flightdata.Ahrs1_bYawRate.data(t_init)];
testttt = lsim(sys,u_input,t);



figure();
subplot(2,2,1);
plot(t, testttt(:,4)+flightdata.Ahrs1_bYawRate.data(t_init), 'Color','g'); hold on;
plot(t, flightdata.Ahrs1_bYawRate.data(t_init:t_end),'Color','b'); hold on;
plot_title = ['Yaw Rate over time - ', num2str(maneuver.Name)];
title(plot_title);
xlabel('Time [s]');
ylabel('Yaw rate [rad/s]');

subplot(2,2,2);
plot(t, testttt(:,3), 'Color','g'); hold on;
plot(t, flightdata.Ahrs1_bRollRate.data(t_init:t_end),'Color','b'); 
plot_title =(['Roll Rate over time - ', num2str(maneuver.Name)]);
title(plot_title);
xlabel('Time [s]');
ylabel('Roll rate [rad/s]');

subplot(2,2,3);
plot(t, testttt(:,2)+flightdata.Ahrs1_Roll.data(t_init), 'Color','g'); hold on;
%plot(t, rad2deg(testttt(:,2)), 'Color','g'); hold on;
plot(t, flightdata.Ahrs1_Roll.data(t_init:t_end),'Color','b');
plot_title = (['Bank Angle over time - ', num2str(maneuver.Name)]);
title(plot_title);
xlabel('Time [s]');
ylabel('Bank Angle [deg]');

subplot(2,2,4);
plot(t, u_input(1,:), 'Color', 'r'); hold on; 
plot(t, u_input(2,:), 'Color', 'm'); 
plot_title = (['Input - ', num2str(maneuver.Name)]);
title(plot_title);
xlabel('Time [s]');
ylabel('Bank Angle [deg]');

filename = ['All_plots_',maneuver.Name,'.png'];
saveas(gcf,filename);

figure('Visible','off');
plot(t, testttt(:,4)+flightdata.Ahrs1_bYawRate.data(t_init), 'Color','g'); hold on;
plot(t, flightdata.Ahrs1_bYawRate.data(t_init:t_end),'Color','b');
plot_title = (['Yaw Rate over time - ', num2str(maneuver.Name)]);
title(plot_title);
xlabel('Time [s]');
ylabel('Yaw rate [rad/s]')
saveas(gcf,[plot_title,'.png']);

figure('Visible','off');
plot(t, testttt(:,3), 'Color','g'); hold on;
plot(t, flightdata.Ahrs1_bRollRate.data(t_init:t_end),'Color','b'); 
plot_title =(['Roll Rate over time - ', num2str(maneuver.Name)]);
title(plot_title);
xlabel('Time [s]');
ylabel('Roll rate [rad/s]');
saveas(gcf,[plot_title,'.png']);

figure('Visible','off');
plot(t, testttt(:,2)+flightdata.Ahrs1_Roll.data(t_init), 'Color','g'); hold on;
plot(t, flightdata.Ahrs1_Roll.data(t_init:t_end),'Color','b');
plot_title = (['Bank Angle over time - ', num2str(maneuver.Name)]);
title(plot_title);
xlabel('Time [s]');
ylabel('Bank Angle [deg]')
saveas(gcf,[plot_title,'.png']);
end