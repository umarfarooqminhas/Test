function [D,R,c,H,Ad,AR,mi,ro,omega] = dataV (N_wg,Nb,HrD,sigma,V0,lam)

D = 2*N_wg-1;     %15  m
R = D/2;          %7.5 m
c = sigma*D/Nb;   %1.25 m
H = HrD * D;      %18  m
AR = H / c ;
Ad = D*H;

mi=1.81*10^(-5);          % Dynamic viscosity           [Kg/ms]
ro=1.22;                  % Rho = density               [kg/m3]
%v=mi / ro;                % Kinematic Viscosity         [m2/s]
%z = 100;                  % altidude	                [m]

rot_rad_s = V0*lam/R;         % Rotational speed        [rad/sec]
%rot_rad_m =	rot_rad_s*60;     % Rotational speed        [rad/min]
%rot_deg_s = rot_rad_s*180/pi; % Rotational speed        [deg/sec]

%rot_rpm = rot_rad_m/(2*pi); % rpm
%rot_rps = rot_rad_s/(2*pi); % rps

omega=rot_rad_s; %it may be changed according to the needs