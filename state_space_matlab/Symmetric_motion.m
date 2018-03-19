% Citation 550 - Linear simulation

load('flighttestdata.mat');
open_data;


%% CALCULATING AIRCRAFT MASS AS A FUNCTION OF TIME
% Calculating the initial mass
m_extra_payload = convmass(220, 'lbm', 'kg');
m_fuel_start = convmass(2650, 'lbm', 'kg');
m_people = 92 + 95 + 76 + 61 + 59 + 66 + 77 + 77 + 84;
m_OEW = convmass(9165, 'lbm', 'kg');

m_start = m_OEW + m_people + m_extra_payload + m_fuel_start;

% Calculate the decrease in mass as a function of time
fuelusedl_kg = convmass (fuelusedl,'lbm', 'kg');
fuelusedr_kg = convmass (fuelusedr, 'lbm', 'kg');

ac_mass_lst=[];

for i =1:length(t)
    
    ac_mass_lst(end+1) = m_start - fuelusedl_kg(i) - fuelusedr_kg(i);
end
%% SYMMETRIC MOTION STATE SPACE INPUTS


% Stationary flight condition

hp0    = h(3200 * 10) * 0.3048;      	  % pressure altitude in the stationary flight condition [m]
V0     = Vtas(3200 * 10) * 0.514444;            % true airspeed in the stationary flight condition [m/sec]
alpha0 = alpha(3200 * 10) * pi/180;       	  % angle of attack in the stationary flight condition [rad]
th0    = pitch(3200 * 10) * pi/180;       	  % pitch angle in the stationary flight condition [rad]

% Aircraft mass
m      = 5500;         	  % mass [kg]

% aerodynamic properties
e      = 0.865;            % Oswald factor [ ]
CD0    = 0.022;            % Zero lift drag coefficient [ ]
CLa    = 4.38;            % Slope of CL-alpha curve [ ]

% Longitudinal stability
Cma    = -0.5626;            % longitudinal stabilty [ ]
Cmde   = -1.1642;            % elevator effectiveness [ ]

% Aircraft geometry

S      = 30.00;	          % wing area [m^2]
Sh     = 0.2*S;           % stabiliser area [m^2]
Sh_S   = Sh/S;	          % [ ]
lh     = 0.71*5.968;      % tail length [m]
c      = 2.0569;	  % mean aerodynamic cord [m]
lh_c   = lh/c;	          % [ ]
b      = 15.911;	  % wing span [m]
bh     = 5.791;	          % stabilser span [m]
A      = b^2/S;           % wing aspect ratio [ ]
Ah     = bh^2/Sh;         % stabilser aspect ratio [ ]
Vh_V   = 1;		  % [ ]
ih     = -2*pi/180;       % stabiliser angle of incidence [rad]

% Constant values concerning atmosphere and gravity

rho0   = 1.2250;          % air density at sea level [kg/m^3] 
lambda = -0.0065;         % temperature gradient in ISA [K/m]
Temp0  = 288.15;          % temperature at sea level in ISA [K]
R      = 287.05;          % specific gas constant [m^2/sec^2K]
g      = 9.81;            % [m/sec^2] (gravity constant)

rho    = rho0*((1+(lambda*hp0/Temp0)))^(-((g/(lambda*R))+1));   % [kg/m^3]  (air density)
W      = m*g;				                        % [N]       (aircraft weight)

% Constant values concerning aircraft inertia

muc    = m/(rho*S*c);
mub    = m/(rho*S*b);
KX2    = 0.019;
KZ2    = 0.042;
KXZ    = 0.002;
KY2    = 1.25*1.114;

% Aerodynamic constants

Cmac   = 0;                     % Moment coefficient about the aerodynamic centre [ ]
CNwa   = CLa;   		        % Wing normal force slope [ ]
CNha   = 2*pi*Ah/(Ah+2);        % Stabiliser normal force slope [ ]
depsda = 4/(A+2);               % Downwash gradient [ ]

% Lift and drag coefficient

CL = 2*W/(rho*V0^2*S);               % Lift coefficient [ ]
CD = CD0 + (CLa*alpha0)^2/(pi*A*e);  % Drag coefficient [ ]

% Stabiblity derivatives

CX0    = W*sin(th0)/(0.5*rho*V0^2*S);
CXu    = -0.02792;
CXa    = -0.47966;
CXadot = +0.08330;
CXq    = -0.28170;
CXde   = -0.03728;

CZ0    = -W*cos(th0)/(0.5*rho*V0^2*S);
CZu    = -0.37616;
CZa    = -5.74340;
CZadot = -0.00350;
CZq    = -5.66290;
CZde   = -0.69612;

Cmu    = +0.06990;
Cmadot = +0.17800;
Cmq    = -8.79415;

CYb    = -0.7500;
CYbdot =  0     ;
CYp    = -0.0304;
CYr    = +0.8495;
CYda   = -0.0400;
CYdr   = +0.2300;

Clb    = -0.10260;
Clp    = -0.71085;
Clr    = +0.23760;
Clda   = -0.23088;
Cldr   = +0.03440;

Cnb    =  +0.1348;
Cnbdot =   0     ;
Cnp    =  -0.0602;
Cnr    =  -0.2061;
Cnda   =  -0.0120;
Cndr   =  -0.0939;

%% Calculate C matrices here

C1_sym = [  -2*muc*c/V0^2        0                               0               0;...
            0                       (CZadot-2*muc)*c/V0    0               0;...
            0                       0                               -c/V0         0;...
            0                       Cmadot*c/V0             0               -2*muc*KY2*c^2/V0^2];
        
        
C2_sym = [CXu     CXa   CZ0       CXq*c/V0;...
            CZu   CZa    -CX0  (CZq+2*muc)*c/V0; ...
            0       0               0       c/V0;...
            Cmu   Cma   0   Cmq*c/V0];
        
C3_sym = [CXde ; CZde ; 0 ; Cmde ];

% Calculate A and B matrices

A = -inv(C1_sym) * C2_sym;
B = -inv(C1_sym) * C3_sym;

% C, D matrices for symmetric motion 

C = eye(5,4);
D = [ 0 ; 0 ; 0 ; 0 ; 1];

%% Create a state space model

sys = ss(A,B,C,D);
time = 0:0.1:60;
u = zeros(size(time));
u(1:10) = delta_e(3200 * 10) * pi/180; 

response = lsim(sys,u,time);

%load('/Users/Mykolas/Documents/GitHub/SVV/Data reading/open.m')
%plot (time, 180/pi*response(:,3),t(0:0.1:310), pitch(3190*10:3500*10) )% t(0:0.1:310), pitch(3190*10:3500*10)) % pitch angle response
%plot ( t(0:0.1:310), pitch(3190*10:3500*10))

