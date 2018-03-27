clear all;
close all;
clc;

load('flightdata.mat')
ac_mass_time;

t_start = 31990;
t_stop = 3960;
t_maneuver = t_start:t_stop;

Cit_par;
V = V0;

da = flightdata.delta_a.data;
dr = flightdata.delta_r.data;

Matrices_a;
u_input = [ (da(t_maneuver)'-da(t_start));...
           -(dr(t_maneuver)'-flightdata.delta_r.data(t_init))];
       
%Clb    = -0.10260;
Clb    = -0.13260;       