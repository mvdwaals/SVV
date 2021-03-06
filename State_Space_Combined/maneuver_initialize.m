%Input 1: Aileron - Input 2: Rudder - Input 3: Elevator

%% Initialize Asymmetric Eigenmotions

%Dutch Roll - Undamped
dutch_roll_undamped = struct('Name','Dutch_Roll_Undamped','t_init',36550,'t_end',36750,'input',2,'mode','a','first_plot','p_roll_rate','second_plot','r_yaw_rate', 'disp', 7);

%Dutch Roll - Yaw Damper
dutch_roll_damped = struct('Name','Dutch_Roll_Damped','t_init',37185,'t_end',37385,'input',2,'mode','a','first_plot','p_roll_rate','second_plot','r_yaw_rate', 'disp', 7);

%Aperiodic Roll
aper_roll = struct('Name','Aperiodic_Roll','t_init',35569,'t_end',35629,'input',1,'mode','a','first_plot','phi_bank_angle','second_plot','p_roll_rate', 'disp', 6);

%Spiral - 1st
spiral_1 = struct('Name','Spiral_1st','t_init',38250,'t_end',38850,'input',2,'mode','a','first_plot','phi_bank_angle','second_plot','r_yaw_rate', 'disp', 9);

%Spiral - 2nd
spiral_2 = struct('Name','Spiral_2nd','t_init',38800,'t_end',39600,'input',2,'mode','a','first_plot','phi_bank_angle','second_plot','r_yaw_rate', 'disp', 9);

maneuvers_asym = [dutch_roll_undamped, dutch_roll_damped, spiral_2, aper_roll];

%% Initialize Symmetric Eigenmotions

%Short Period
short_period = struct('Name','Short_Period','t_init',31990,'t_end',32015,'input',3,'mode','s','first_plot','theta_pitch_angle','second_plot','q_pitch_rate', 'disp', 8);

%Phugoid
phugoid = struct('Name','Phugoid','t_init',32477,'t_end',32677,'input',3,'mode','s','first_plot','theta_pitch_angle','second_plot','V', 'disp', 5);

maneuvers_sym = [short_period, phugoid];

%Group all maneuvers

maneuvers = [dutch_roll_undamped, dutch_roll_damped, spiral_2, aper_roll, short_period, phugoid];