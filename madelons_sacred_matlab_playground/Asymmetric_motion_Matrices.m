

Asymmetric_motion;
V=V0;

C1_asym = [(CYbdot-2*mub)*b/V    	0       0       0;...
            0                   -1/2*b/V	0	0;...
            0                   0           -4*mub*KX2*b^2/(2*V^2)	-4*mub*KXZ*b^2/(2*V^2);...
            Cnbdot*b/V          0           -4*mub*KXZ*b^2/(2*V^2)  -4*mub*KZ2*b^2/(2*V^2)];
        
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
 
t = flightdata.time.data(35500:36000);
 
 u_input = [flightdata.delta_a.data(35500:36000)'-flightdata.delta_a.data(35500);...
            flightdata.delta_r.data(35500:36000)'-flightdata.delta_r.data(35500)];
     
 
 sys = ss(A,B,C,D);

 
 
testttt = lsim(sys,u_input,t);

figure(1);
plot(t, testttt(:,4)+flightdata.Ahrs1_bYawRate.data(35500), t, flightdata.Ahrs1_bYawRate.data(35500:36000)); hold off;
title('Yaw Rate over time - Aperiodic Roll');
xlabel('');
ylabel('');

figure(2);
plot(t, testttt(:,3), t, flightdata.Ahrs1_bRollRate.data(35500:36000));
title('Roll Rate over time - Aperiodic Roll');
xlabel('');
ylabel('');