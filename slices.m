function [N_st,D_theta,theta_v_rad_up,theta_v_rad_down] = slices(N_slices)

N_st = N_slices/2;
D_theta_deg = 360/N_slices;
D_theta_rad = D_theta_deg*pi/180;
D_theta = D_theta_rad;

%N_st_v = 1:1:N_st_tot;

%upstream

theta_in_up = 0;
theta_end_up = 180;
%theta_v_deg = theta_in:D_theta_deg:theta_end;
theta_v_deg_up = theta_in_up+D_theta_deg/2:D_theta_deg:theta_end_up-D_theta_deg/2;

theta_v_deg_up = theta_v_deg_up';
theta_v_rad_up = theta_v_deg_up.*pi./180;


%downstream

theta_in_down = 360;
theta_end_down = 180;
%theta_v_deg = theta_in:D_theta_deg:theta_end;
theta_v_deg_down = theta_in_down-D_theta_deg/2:-D_theta_deg:theta_end_down+D_theta_deg/2;

theta_v_deg_down = theta_v_deg_down';
theta_v_rad_down = theta_v_deg_down.*pi./180;
