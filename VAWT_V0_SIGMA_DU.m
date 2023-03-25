%% VAWT AT FIXED V0 = 10 m/s
clear 
close all
clc

set(0,'DefaultAxesFontsize',10);                                           % grandezza numeri assi
set(0,'Defaultlinelinewidth',1.5);                                           % spessore linea tracciato

%% LOADING and EXTRACTION

[NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M,DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M] = loading ();

%% SLICES

N_slices = 40;
[N_st,D_theta,theta_v_rad_up,theta_v_rad_down] = slices(N_slices);

%% CX VECTOR

Cx_a = 1.816; %Cx_a = 3.5; 
Cx_a = 2;
at = 1 - 0.5 * Cx_a^(0.5);
step_a = 0.001;

[a_v, Cx_v] = Cx_vector (Cx_a,at,step_a);

figure(1)
plot(a_v,Cx_v,'b')
xlabel('a')
ylabel('Cx(a)')

N_wg = 8;
Nb = 3;
HrD = 1.2;

V0 = 10;                   % Wind velocity (undisturbed) [m/s]

s0 = 0.250;
s1 = 0.125;
s2 = 0.375;
s3 = 0.500;

%% 

%% LAMBDA = 2

sigma = s1;
lam = 2;                  % lambda = tip speed ratio    [-] %5
[D,R,c,H,Ad,AR,mi,ro,omega] = dataV (N_wg,Nb,HrD,sigma,V0,lam); % omega=rot_rad/s

[a1_v_L2,T_D1,TAO_1] = VAWT_UP_DU (N_st, D_theta, theta_v_rad_up, a_v, Cx_v, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);
[a2_v_L2_s1,T_D2,TAO_2] = VAWT_DOWN_DU_pp_0 (N_st, D_theta, theta_v_rad_down,a1_v_L2, at, Cx_a, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);

P_id_L2 = 0.5 * ro * V0^3 * Ad ;
P_av_L2 = (Nb/(2*Ad)) * ( sum(TAO_1) + sum(TAO_2));
Cp_L2_s1 = P_av_L2 / P_id_L2

%% lambda = 2 & sigma 0

sigma = s0;
lam = 2;                  % lambda = tip speed ratio    [-] %5

[D,R,c,H,Ad,AR,mi,ro,omega] = dataV (N_wg,Nb,HrD,sigma,V0,lam); % omega=rot_rad/s

[a1_v_L2,T_D1,TAO_1] = VAWT_UP_DU (N_st, D_theta, theta_v_rad_up, a_v, Cx_v, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);
[a2_v_L2_s0,T_D2,TAO_2] = VAWT_DOWN_DU_pp_0 (N_st, D_theta, theta_v_rad_down,a1_v_L2, at, Cx_a, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);

%% Cp for each sigma at lambda = 2

sigma_v = s0:0.025:s3;

a2_v_L2_0 = a2_v_L2_s0;

for i = 1 : length(sigma_v)

    sigma = sigma_v(i);
    [D,R,c,H,Ad,AR,mi,ro,omega] = dataV (N_wg,Nb,HrD,sigma,V0,lam);

    [a1_v_L2,T_D1,TAO_1] = VAWT_UP_DU (N_st, D_theta, theta_v_rad_up, a_v, Cx_v, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);
    [a2_v_L2_s2,T_D2,TAO_2] = VAWT_DOWN_DU_pp (N_st, D_theta, theta_v_rad_down, a1_v_L2,a2_v_L2_0, at, Cx_a, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);

    P_id_L2 = 0.5 * ro * V0^3 * Ad ;
    P_av_L2 = (Nb/(2*Ad)) * ( sum(TAO_1) + sum(TAO_2));
    Cp_L2_v(i) = P_av_L2 / P_id_L2;

    a2_v_L2_0 = a2_v_L2_s2;

end

sigma_v = [s1,sigma_v];
Cp_L2_v = [Cp_L2_s1,Cp_L2_v];

figure(2)
plot(sigma_v,Cp_L2_v)
ylabel('Cp at lambda = 2')
xlabel('solidity')

%%

%% LAMBDA = 3

Cx_a = 1.816; %Cx_a = 3.5; 
Cx_a = 2;
at = 1 - 0.5 * Cx_a^(0.5);
[a_v, Cx_v] = Cx_vector (Cx_a,at,step_a);

sigma = s1;
lam = 3;                  % lambda = tip speed ratio    [-] %5
[D,R,c,H,Ad,AR,mi,ro,omega] = dataV (N_wg,Nb,HrD,sigma,V0,lam); % omega=rot_rad/s

[a1_v_L3,T_D1,TAO_1] = VAWT_UP_DU (N_st, D_theta, theta_v_rad_up, a_v, Cx_v, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);
[a2_v_L3_s1,T_D2,TAO_2] = VAWT_DOWN_DU_pp_0 (N_st, D_theta, theta_v_rad_down,a1_v_L3, at, Cx_a, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);
%[a2_v_L,T_D2,TAO_2] = VAWT_DOWN_DU_pp (N_st, D_theta, theta_v_rad_down, a1_v_L4,a2_v_L3, at, Cx_a_4, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);

P_id_L3 = 0.5 * ro * V0^3 * Ad ;
P_av_L3 = (Nb/(2*Ad)) * ( sum(TAO_1) + sum(TAO_2));
Cp_L3_s1 = P_av_L3 / P_id_L3


sigma = s0;
lam = 3;                  % lambda = tip speed ratio    [-] %5
[D,R,c,H,Ad,AR,mi,ro,omega] = dataV (N_wg,Nb,HrD,sigma,V0,lam); % omega=rot_rad/s

[a1_v_L3,T_D1,TAO_1] = VAWT_UP_DU (N_st, D_theta, theta_v_rad_up, a_v, Cx_v, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);
[a2_v_L3_s0,T_D2,TAO_2] = VAWT_DOWN_DU_pp_0 (N_st, D_theta, theta_v_rad_down,a1_v_L3, at, Cx_a, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);

%%  Cp for each sigma at lambda = 3

sigma_v = s0:0.025:s3;

a2_v_L3_0 = a2_v_L3_s0;

for i = 1 : length(sigma_v)

    sigma = sigma_v(i);
    [D,R,c,H,Ad,AR,mi,ro,omega] = dataV (N_wg,Nb,HrD,sigma,V0,lam);

    [a1_v_L3,T_D1,TAO_1] = VAWT_UP_DU (N_st, D_theta, theta_v_rad_up, a_v, Cx_v, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);
    [a2_v_L3_s2,T_D2,TAO_2] = VAWT_DOWN_DU_pp (N_st, D_theta, theta_v_rad_down, a1_v_L3,a2_v_L3_0, at, Cx_a, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);

    P_id_L3 = 0.5 * ro * V0^3 * Ad ;
    P_av_L3 = (Nb/(2*Ad)) * ( sum(TAO_1) + sum(TAO_2));
    Cp_L3_v(i) = P_av_L3 / P_id_L3;

    a2_v_L3_0 = a2_v_L3_s2;

end

sigma_v = [s1,sigma_v];
Cp_L3_v = [Cp_L3_s1,Cp_L3_v];

figure(3)
plot(sigma_v,Cp_L3_v)
ylabel('Cp at lambda = 3')
xlabel('solidity')

%%

%% LAMBDA = 4

Cx_a = 2;
at = 1 - 0.5 * Cx_a^(0.5);
[a_v, Cx_v] = Cx_vector (Cx_a,at,step_a);


sigma = s1;                  % Wind velocity (undisturbed) [m/s] % 7.3
lam = 4;                  % lambda = tip speed ratio    [-] %5
[D,R,c,H,Ad,AR,mi,ro,omega] = dataV (N_wg,Nb,HrD,sigma,V0,lam); % omega=rot_rad/s

[a1_v_L4,T_D1,TAO_1] = VAWT_UP_DU (N_st, D_theta, theta_v_rad_up, a_v, Cx_v, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);
[a2_v_L4_s1,T_D2,TAO_2] = VAWT_DOWN_DU_pp_0 (N_st, D_theta, theta_v_rad_down,a1_v_L4, at, Cx_a, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);
%[a2_v_L,T_D2,TAO_2] = VAWT_DOWN_DU_pp (N_st, D_theta, theta_v_rad_down, a1_v_L4,a2_v_L3, at, Cx_a_4, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);
P_id_L4 = 0.5 * ro * V0^3 * Ad ;
P_av_L4 = (Nb/(2*Ad)) * ( sum(TAO_1) + sum(TAO_2));
Cp_L4_s1 = P_av_L4 / P_id_L4


% SIGMA 0 % lambda 3
sigma = s0;                  % Wind velocity (undisturbed) [m/s] % 7.3
lam = 3;                  % lambda = tip speed ratio    [-] %5
[D,R,c,H,Ad,AR,mi,ro,omega] = dataV (N_wg,Nb,HrD,sigma,V0,lam); % omega=rot_rad/s

[a1_v_L3,T_D1,TAO_1] = VAWT_UP_DU (N_st, D_theta, theta_v_rad_up, a_v, Cx_v, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);
[a2_v_L3_s0,T_D2,TAO_2] = VAWT_DOWN_DU_pp_0 (N_st, D_theta, theta_v_rad_down,a1_v_L3, at, Cx_a, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);
          
% SIGMA 0 % lambda 4
Cx_a_4 = 2.7; %Cx_a = 1.816;
at_4 = 1 - 0.5 * Cx_a_4^(0.5);
[a_v_4, Cx_v_4] = Cx_vector (Cx_a_4,at_4,step_a); %figure plot(a_v,Cx_v,'b') xlabel('a') ylabel('Cx(a)')

sigma = s0; 
lam = 4;                  % lambda = tip speed ratio    [-] %5
[D,R,c,H,Ad,AR,mi,ro,omega] = dataV (N_wg,Nb,HrD,sigma,V0,lam); % omega=rot_rad/s

[a1_v_L4,T_D1,TAO_1] = VAWT_UP_DU (N_st, D_theta, theta_v_rad_up, a_v_4, Cx_v_4, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);
[a2_v_L4_s0,T_D2,TAO_2]  = VAWT_DOWN_DU_pp (N_st, D_theta, theta_v_rad_down, a1_v_L4,a2_v_L3_s0, at_4, Cx_a_4, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);


%% Cp for each sigma at lambda = 4

sigma_v = s0:0.025:s3;

a2_v_L4_0 = a2_v_L4_s0;

for i = 1 : length(sigma_v)

    sigma = sigma_v(i);
    [D,R,c,H,Ad,AR,mi,ro,omega] = dataV (N_wg,Nb,HrD,sigma,V0,lam);

    [a1_v_L4,T_D1,TAO_1] = VAWT_UP_DU (N_st, D_theta, theta_v_rad_up, a_v, Cx_v, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);
    [a2_v_L4_s2,T_D2,TAO_2] = VAWT_DOWN_DU_pp (N_st, D_theta, theta_v_rad_down, a1_v_L4,a2_v_L4_0, at, Cx_a, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);

    P_id_L4 = 0.5 * ro * V0^3 * Ad ;
    P_av_L4 = (Nb/(2*Ad)) * ( sum(TAO_1) + sum(TAO_2));
    Cp_L4_v(i) = P_av_L4 / P_id_L4;

    a2_v_L4_0 = a2_v_L4_s2;

end

sigma_v = [s1,sigma_v];
Cp_L4_v = [Cp_L4_s1,Cp_L4_v];

figure(4)
plot(sigma_v,Cp_L4_v)
ylabel('Cp at lambda = 4')
xlabel('solidity')

%%

%% lambda = 5 

Cx_a = 1.816; %Cx_a = 3.5;
Cx_a = 2;
at = 1 - 0.5 * Cx_a^(0.5);
[a_v, Cx_v] = Cx_vector (Cx_a,at,step_a);

sigma = s1;               
lam = 5;                  % lambda = tip speed ratio    [-] %5
[D,R,c,H,Ad,AR,mi,ro,omega] = dataV (N_wg,Nb,HrD,sigma,V0,lam); % omega=rot_rad/s

[a1_v_L5,T_D1,TAO_1] = VAWT_UP_DU (N_st, D_theta, theta_v_rad_up, a_v, Cx_v, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);
[a2_v_L5_s1,T_D2,TAO_2] = VAWT_DOWN_DU_pp_0 (N_st, D_theta, theta_v_rad_down,a1_v_L5, at, Cx_a, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);
%[a2_v_L,T_D2,TAO_2] = VAWT_DOWN_DU_pp (N_st, D_theta, theta_v_rad_down, a1_v_L4,a2_v_L3, at, Cx_a_4, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);

P_id_L5 = 0.5 * ro * V0^3 * Ad ;
P_av_L5 = (Nb/(2*Ad)) * ( sum(TAO_1) + sum(TAO_2));
Cp_L5_s1 = P_av_L5 / P_id_L5


% lambda = 3
sigma = s0;                 
lam = 3;                  % lambda = tip speed ratio    [-] %5
[D,R,c,H,Ad,AR,mi,ro,omega] = dataV (N_wg,Nb,HrD,sigma,V0,lam); % omega=rot_rad/s

[a1_v_L3,T_D1,TAO_1] = VAWT_UP_DU (N_st, D_theta, theta_v_rad_up, a_v, Cx_v, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);
[a2_v_L3_s0,T_D2,TAO_2] = VAWT_DOWN_DU_pp_0 (N_st, D_theta, theta_v_rad_down,a1_v_L3, at, Cx_a, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);
    

% lambda = 4
Cx_a_4 = 2.7; %Cx_a = 1.816;
at_4 = 1 - 0.5 * Cx_a_4^(0.5);
[a_v_4, Cx_v_4] = Cx_vector (Cx_a_4,at_4,step_a); %figure plot(a_v,Cx_v,'b') xlabel('a') ylabel('Cx(a)')

sigma = s0;
lam = 4;                  % lambda = tip speed ratio    [-] %5
[D,R,c,H,Ad,AR,mi,ro,omega] = dataV (N_wg,Nb,HrD,sigma,V0,lam); % omega=rot_rad/s

[a1_v_L4,T_D1,TAO_1] = VAWT_UP_DU (N_st, D_theta, theta_v_rad_up, a_v, Cx_v, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);
[a2_v_L4_s0,T_D2,TAO_2]  = VAWT_DOWN_DU_pp (N_st, D_theta, theta_v_rad_down, a1_v_L4,a2_v_L3_s0, at_4, Cx_a_4, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);


% lambda = 5
Cx_a_5 = 3.3; %Cx_a = 1.816;
at_5 = 1 - 0.5 * Cx_a_5^(0.5);
[a_v_5, Cx_v_5] = Cx_vector (Cx_a_5,at_5,step_a); %figure plot(a_v,Cx_v,'b') xlabel('a') ylabel('Cx(a)')

sigma = s0;
lam = 5;                  % lambda = tip speed ratio    [-] %5
[D,R,c,H,Ad,AR,mi,ro,omega] = dataV (N_wg,Nb,HrD,sigma,V0,lam); % omega=rot_rad/s

[a1_v_L5,T_D1,TAO_1] = VAWT_UP_DU (N_st, D_theta, theta_v_rad_up, a_v_5, Cx_v_5, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);
[a2_v_L5_s0,T_D2,TAO_2] = VAWT_DOWN_DU_pp (N_st, D_theta, theta_v_rad_down, a1_v_L5,a2_v_L4_s0, at_5, Cx_a_5, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);

%% Cp for each sigma at lambda = 5

sigma_v = s0:0.025:s3;

a2_v_L5_0 = a2_v_L5_s0;

for i = 1 : length(sigma_v)

    sigma = sigma_v(i);
    [D,R,c,H,Ad,AR,mi,ro,omega] = dataV (N_wg,Nb,HrD,sigma,V0,lam);

    [a1_v_L5,T_D1,TAO_1] = VAWT_UP_DU (N_st, D_theta, theta_v_rad_up, a_v, Cx_v, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);
    [a2_v_L5_s2,T_D2,TAO_2] = VAWT_DOWN_DU_pp (N_st, D_theta, theta_v_rad_down, a1_v_L5,a2_v_L5_0, at, Cx_a, Nb, c,ro,mi,H,R,AR, V0, lam, omega, DU_40K,DU_80K,DU_160K,DU_360K,DU_700K,DU_1M,DU_2M,DU_5M);

    P_id_L5 = 0.5 * ro * V0^3 * Ad ;
    P_av_L5 = (Nb/(2*Ad)) * ( sum(TAO_1) + sum(TAO_2));
    Cp_L5_v(i) = P_av_L5 / P_id_L5;

    a2_v_L5_0 = a2_v_L5_s2;

end

sigma_v = [s1,sigma_v];
Cp_L5_v = [Cp_L5_s1,Cp_L5_v];

figure(5)
plot(sigma_v,Cp_L5_v)
ylabel('Cp at lambda = 5')
xlabel('solidity')

%%

%% images

sigma_v = s0:0.025:s3;
sigma_v = [s1, sigma_v];

lambda_v = 2:1:5;

jj = [1 2 7 ];%12];

for i = 1:length(jj)
    z = jj(i);
    j2 = Cp_L2_v(z);
    j3 = Cp_L3_v(z);
    j4 = Cp_L4_v(z);
    j5 = Cp_L5_v(z);

    Cp_s_v = [j2 j3 j4 j5];

    figure(7)
    plot(lambda_v,Cp_s_v,'-o')
    hold on
    xlabel('lambda [-] ')
    ylabel('Cp [-] ')
    ylim([0 0.7])
    legend('sigma = 0.125','sigma = 0.250', 'sigma = 0.375', 'sigma = 0.500','Location','southwest')


end

%% images

Cp_L2_vec = [0.271141515660505	0.454505666120728	0.456678673880480	0.446684348033421	0.430317645346600	0.416062064130531	0.401253549031664	0.394473845360590	0.379271684636639	0.355076506281002	0.334053318195181	0.313306718408368];

Cp_L3_vec = [0.612420897311369	0.479282002828426	0.448947635273489	0.419289383326765	0.391161416402762	0.362923648237405	0.334838862751294	0.308038004943058	0.281439397179873	0.255070921182718	0.229861883417015	0.206444119098801];
  
Cp_L4_vec = [0.591370525395186	0.389358360656218	0.350254177738281	0.312627059842922	0.278275733171201	0.242674532487097	0.211859215957643	0.181435199963635	0.149004512063094	0.120579055812827	0.0958294053842310	0.0652063156082660];

Cp_L5_vec = [0.264813256254312	0.222802530111233	0.188764605380991	0.160817580975115	0.120152978441176	0.0889925629171390	0.0562422720404380	0.0154274477222010 0 0 0 0];
    %0.554795621623598 0.334306586552559	0.281782682992607	0.235189345917949	0.191023114307577	0.148653388258926	0.108177974512491	0.0760052698195887	0.0409103256753319	0.00138831193398600	-0.0360654957412800	-0.0719803177373988];

%{
figure plot(Cp_L2_vec, '-o') hold on plot(Cp_L3_vec, '-o') hold on plot(Cp_L4_vec, '-o') hold on plot(Cp_L5_vec, '-o')
%}

sigma_v = s0:0.025:s3;
sigma_v = [s1, sigma_v];

lambda_v = 2:1:5;

jj = [1 2 7];% 12];

for i = 1:length(jj)
    z = jj(i);
    j2 = Cp_L2_vec(z);
    j3 = Cp_L3_vec(z);
    j4 = Cp_L4_vec(z);
    j5 = Cp_L5_vec(z);

    Cp_s_v_saved = [j2 j3 j4 j5];

    figure(8)
    plot(lambda_v,Cp_s_v_saved,'-o')
    hold on
    xlabel('lambda [-] ')
    ylabel('Cp [-] ')
    ylim([0 0.7])
    legend('sigma = 0.125','sigma = 0.250', 'sigma = 0.375', 'sigma = 0.500','Location','southwest')

end


lambda_v = 2:1:4;
jj = 12;

for i = 1:length(jj)
    z = jj(i);
    j2 = Cp_L2_vec(z);
    j3 = Cp_L3_vec(z);
    j4 = Cp_L4_vec(z);
    %j5 = Cp_L5_vec(z);

    Cp_s_v_saved = [j2 j3 j4 ];

    figure(8)
    plot(lambda_v,Cp_s_v_saved,'-o')
    hold on
    xlabel('lambda [-] ')
    ylabel('Cp [-] ')
    ylim([0 0.7])
    legend('sigma = 0.125','sigma = 0.250', 'sigma = 0.375', 'sigma = 0.500','Location','southwest')

end



