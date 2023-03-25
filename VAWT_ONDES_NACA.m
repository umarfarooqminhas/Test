%% VAWT ON DESIGN
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
%Cx_a = 2;
at = 1 - 0.5 * Cx_a^(0.5);
step_a = 0.001;

[a_v, Cx_v] = Cx_vector (Cx_a,at,step_a);

figure(1)
plot(a_v,Cx_v,'b')
xlabel('a')
ylabel('Cx(a)')

%% CONSTRAINTS

N_wg = 8;
Nb = 3;
HrD = 1.2;

%% lambda = 2

sigma = 0.25;
V0 = 5;                   % Wind velocity (undisturbed) [m/s] % 7.3
lam = 2;                  % lambda = tip speed ratio    [-] %5

[D,R,c,H,Ad,AR,mi,ro,omega] = dataV (N_wg,Nb,HrD,sigma,V0,lam); % omega=rot_rad/s

[a1_v_L2,T_D1,TAO_1] = VAWT_UP_NACA (N_st, D_theta, theta_v_rad_up, a_v, Cx_v, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M);
[a2_v_L2,T_D2,TAO_2] = VAWT_DOWN_NACA_pp_0 (N_st, D_theta, theta_v_rad_down,a1_v_L2, at, Cx_a, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M);

P_id_L2 = 0.5 * ro * V0^3 * Ad ;
P_av_L2 = (Nb/(2*Ad)) * ( sum(TAO_1) + sum(TAO_2));
Cp_L2 = P_av_L2 / P_id_L2;

a1_v = a1_v_L2;
a2_v = a2_v_L2;

figure(2)
plot(a1_v,'r')
hold on
plot(a2_v,'m')
grid on
title('a_D_1 & a_D_2 for each stream tube')
xlabel('N of stream tube')
ylabel('a_D')
legend('a_D_1', 'a_D_2', 'location', 'Northeast')


%% lambda = 3

sigma = 0.25;
V0 = 5;                   % Wind velocity (undisturbed) [m/s] % 7.3
lam = 3;                  % lambda = tip speed ratio    [-] %5

[D,R,c,H,Ad,AR,mi,ro,omega] = dataV (N_wg,Nb,HrD,sigma,V0,lam); % omega=rot_rad_s

[a1_v_L3,T_D1,TAO_1] = VAWT_UP_NACA (N_st, D_theta, theta_v_rad_up, a_v, Cx_v, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M);
[a2_v_L3,T_D2,TAO_2] = VAWT_DOWN_NACA_pp_0 (N_st, D_theta, theta_v_rad_down,a1_v_L3, at, Cx_a, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M);

P_id_L3 = 0.5 * ro * V0^3 * Ad ;
P_av_L3 = (Nb/(2*Ad)) * ( sum(TAO_1) + sum(TAO_2));
Cp_L3 = P_av_L3 / P_id_L3;

a1_v = a1_v_L3;
a2_v = a2_v_L3;

figure(3)
plot(a1_v,'r')
hold on
plot(a2_v,'m')
grid on
title('a_D_1 & a_D_2 for each stream tube')
xlabel('N of stream tube')
ylabel('a_D')
legend('a_D_1', 'a_D_2', 'location', 'Northeast')

%% lambda = 4

Cx_a_4 = 2.8; %Cx_a = 1.816;
at_4 = 1 - 0.5 * Cx_a_4^(0.5);
step_a = 0.001;

[a_v_4, Cx_v_4] = Cx_vector (Cx_a_4,at_4,step_a); %figure plot(a_v,Cx_v,'b') xlabel('a') ylabel('Cx(a)')

sigma = 0.25;
V0 = 5;                   % Wind velocity (undisturbed) [m/s] % 7.3
lam = 4;                  % lambda = tip speed ratio    [-] %5

[D,R,c,H,Ad,AR,mi,ro,omega] = dataV (N_wg,Nb,HrD,sigma,V0,lam); % omega=rot_rad_s

%[a1_v_L4,T_D1,TAO_1] = VAWT_UP_NACA (N_st, D_theta, theta_v_rad_up, a_v, Cx_v, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M);
[a1_v_L4,T_D1,TAO_1] = VAWT_UP_NACA (N_st, D_theta, theta_v_rad_up, a_v_4, Cx_v_4, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M);
[a2_v_L4,T_D2,TAO_2] = VAWT_DOWN_NACA_pp (N_st, D_theta, theta_v_rad_down, a1_v_L4,a2_v_L3, at_4, Cx_a_4, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M);

P_id_L4 = 0.5 * ro * V0^3 * Ad ;
P_av_L4 = (Nb/(2*Ad)) * ( sum(TAO_1) + sum(TAO_2));
Cp_L4 = P_av_L4 / P_id_L4;

a1_v = a1_v_L4;
a2_v = a2_v_L4;

figure(4)
plot(a1_v,'r')
hold on
plot(a2_v,'m')
grid on
title('a_D_1 & a_D_2 for each stream tube')
xlabel('N of stream tube')
ylabel('a_D')
legend('a_D_1', 'a_D_2', 'location', 'Northeast')

%% lambda = 5

Cx_a_5 = 3.3; %Cx_a = 1.816;
at_5 = 1 - 0.5 * Cx_a_5^(0.5);
step_a = 0.001;

[a_v_5, Cx_v_5] = Cx_vector (Cx_a_5,at_5,step_a); %figure plot(a_v,Cx_v,'b') xlabel('a') ylabel('Cx(a)')

sigma = 0.25;
V0 = 5;                  % Wind velocity (undisturbed) [m/s] % 7.3
lam = 5;                  % lambda = tip speed ratio    [-] %5

[D,R,c,H,Ad,AR,mi,ro,omega] = dataV (N_wg,Nb,HrD,sigma,V0,lam); % omega=rot_rad_s

%[a1_v_L5,T_D1,TAO_1] = VAWT_UP_NACA (N_st, D_theta, theta_v_rad_up, a_v, Cx_v, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M);
[a1_v_L5,T_D1,TAO_1] = VAWT_UP_NACA (N_st, D_theta, theta_v_rad_up, a_v_5, Cx_v_5, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M);
[a2_v_L5,T_D2,TAO_2] = VAWT_DOWN_NACA_pp (N_st, D_theta, theta_v_rad_down, a1_v_L5,a2_v_L4, at_5, Cx_a_5, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M);

P_id_L5 = 0.5 * ro * V0^3 * Ad ;
P_av_L5 = (Nb/(2*Ad)) * ( sum(TAO_1) + sum(TAO_2));
Cp_L5 = P_av_L5 / P_id_L5;

a1_v = a1_v_L5;
a2_v = a2_v_L5;

figure(5)
plot(a1_v,'r')
hold on
plot(a2_v,'m')
grid on
title('a_D_1 & a_D_2 for each stream tube')
xlabel('N of stream tube')
ylabel('a_D')
legend('a_D_1', 'a_D_2', 'location', 'Northeast') 

%% FIGURES

%a1_v = a1_v_L5; a2_v = a2_v_L5;

%figure plot(a1_v,'r')
%hold on plot(a2_v,'m') grid on title('a_D_1 & a_D_2 for each stream tube')
%xlabel('N of stream tube') ylabel('a_D') legend('a_D_1', 'a_D_2', 'location', 'Northeast')
%fprintf('a_D_1 and a_D_2 for V0 = %f and lam = %f \n', V0, lam)

%% POWER AND cp
%P_id = 0.5 * ro * V0^3 * Ad ; P_av = (Nb/(2*Ad)) * ( sum(TAO_1) + sum(TAO_2)); Cp = P_av / P_id;

%% Cx calculation

a1_v = a1_v_L3;
a2_v = a2_v_L3;

for i = 1 : N_st
    Cx_1(i)= 4*a1_v(i)*(1-a1_v(i));
    Cx_2(i)= 4*a2_v(i)*(1-a2_v(i))*(1-a1_v(i))^2;
    Cx_tot(i) = Cx_1(i) + Cx_2(i);
end

figure
plot(Cx_1, '-o')
hold on
plot(Cx_2, '-o')
hold on
plot(Cx_tot)
xlabel('N of stream tube')
ylabel('Cx')

%%

a_v_l2 = [a1_v_L2, fliplr(a2_v_L2)];
a_v_l3 = [a1_v_L3, fliplr(a2_v_L3)];
a_v_l4 = [a1_v_L4, fliplr(a2_v_L4)];
a_v_l5 = [a1_v_L5, fliplr(a2_v_L5)];

theta_v_rad_up = theta_v_rad_up';
theta_v_rad_down = theta_v_rad_down';

theta_v=[theta_v_rad_up, fliplr(theta_v_rad_down)];

figure
plot(theta_v,a_v_l2')
hold on
plot(theta_v,a_v_l3')
hold on
plot(theta_v,a_v_l4')
hold on
plot(theta_v,a_v_l5')
xlabel('N of stream tube')
ylabel('a_D')


