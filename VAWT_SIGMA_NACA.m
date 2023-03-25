%% VAWT SIGMA FIXED
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

Cx_a = 1.816; 
%Cx_a = 3.5;
at = 1 - 0.5 * Cx_a^(0.5);
step_a = 0.001;

[a_v, Cx_v] = Cx_vector (Cx_a,at,step_a);

%figure plot(a_v,Cx_v,'b') xlabel('a') ylabel('Cx(a)')

%% DATA

N_wg = 8;
Nb = 3;
HrD = 1.2;

sigma = 0.25;
V0 = 5;                   % Wind velocity (undisturbed) [m/s] % 7.3
lam = 2;                  % lambda = tip speed ratio    [-] %5
[D,R,c,H,Ad,AR,mi,ro,omega] = dataV (N_wg,Nb,HrD,sigma,V0,lam); % omega=rot_rad/s

[a1_v_L2,T_D1,TAO_1] = VAWT_UP_NACA (N_st, D_theta, theta_v_rad_up, a_v, Cx_v, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M);
[a2_v_L2,T_D2,TAO_2] = VAWT_DOWN_NACA_pp_0 (N_st, D_theta, theta_v_rad_down,a1_v_L2, at, Cx_a, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M);

V0 = 5;                   % Wind velocity (undisturbed) [m/s] % 7.3
lam = 3;   
[D,R,c,H,Ad,AR,mi,ro,omega] = dataV (N_wg,Nb,HrD,sigma,V0,lam);
[a1_v_L3,T_D1,TAO_1] = VAWT_UP_NACA (N_st, D_theta, theta_v_rad_up, a_v, Cx_v, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M);
[a2_v_L3,T_D2,TAO_2] = VAWT_DOWN_NACA_pp_0 (N_st, D_theta, theta_v_rad_down,a1_v_L3, at, Cx_a, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M);

%%

Cx_a_4 = 2.7; %Cx_a = 1.816;
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
%a2_v_L4=a2_v_L4.*0.85;

%%

Cx_a_5 = 3.3; %Cx_a = 1.816;
at_5 = 1 - 0.5 * Cx_a_5^(0.5);
step_a = 0.001;
[a_v_5, Cx_v_5] = Cx_vector (Cx_a_5,at_5,step_a); %figure plot(a_v,Cx_v,'b') xlabel('a') ylabel('Cx(a)')
sigma = 0.25;
V0 = 5;                   % Wind velocity (undisturbed) [m/s] % 7.3
lam = 5;                  % lambda = tip speed ratio    [-] %5
[D,R,c,H,Ad,AR,mi,ro,omega] = dataV (N_wg,Nb,HrD,sigma,V0,lam); % omega=rot_rad_s
%[a1_v_L5,T_D1,TAO_1] = VAWT_UP_NACA (N_st, D_theta, theta_v_rad_up, a_v, Cx_v, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M);
[a1_v_L5,T_D1,TAO_1] = VAWT_UP_NACA (N_st, D_theta, theta_v_rad_up, a_v_5, Cx_v_5, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M);
[a2_v_L5,T_D2,TAO_2] = VAWT_DOWN_NACA_pp (N_st, D_theta, theta_v_rad_down, a1_v_L5,a2_v_L4, at_5, Cx_a_5, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M);
%a2_v_L5=a2_v_L5.*0.85;

%%

sigma = 0.25;

V0_min = 5;
V0_Max = 15;

lam_min = 2;
lam_Max = 5;

V0_v = V0_min: 2: V0_Max;
lam_v = lam_min: 1: lam_Max;

for jv = 1 : length(V0_v)
    for jl = 1 : length(lam_v)

        V0 = V0_v(jv);
        lam = lam_v(jl);
        [D,R,c,H,Ad,AR,mi,ro,omega] = dataV (N_wg,Nb,HrD,sigma,V0,lam);

        if lam == 2
            [a1_v,T_D1,TAO_1] = VAWT_UP_NACA (N_st, D_theta, theta_v_rad_up, a_v, Cx_v, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M);
            [a2_v,T_D2,TAO_2] = VAWT_DOWN_NACA_pp_0 (N_st, D_theta, theta_v_rad_down,a1_v, at, Cx_a, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M);

        elseif lam == 3
            [a1_v,T_D1,TAO_1] = VAWT_UP_NACA (N_st, D_theta, theta_v_rad_up, a_v, Cx_v, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M);
            [a2_v,T_D2,TAO_2] = VAWT_DOWN_NACA_pp_0 (N_st, D_theta, theta_v_rad_down,a1_v, at, Cx_a, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M);

        elseif lam == 4
            [a1_v,T_D1,TAO_1] = VAWT_UP_NACA (N_st, D_theta, theta_v_rad_up, a_v_4, Cx_v_4, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M);
            [a2_v,T_D2,TAO_2] = VAWT_DOWN_NACA_pp (N_st, D_theta, theta_v_rad_down, a1_v,a2_v_L3, at_4, Cx_a_4, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M);

        else
            [a1_v,T_D1,TAO_1] = VAWT_UP_NACA (N_st, D_theta, theta_v_rad_up, a_v_5, Cx_v_5, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M);
            [a2_v,T_D2,TAO_2] = VAWT_DOWN_NACA_pp (N_st, D_theta, theta_v_rad_down, a1_v,a2_v_L4, at_5, Cx_a_5, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M);

        end

        P_id = 0.5 * ro * V0^3 * Ad ;
        P_id = P_id / 1000;
        P_av = (Nb/(2*Ad)) * ( sum(TAO_1) + sum(TAO_2));
        P_v(jl) = P_av/1000;
        Cp_v(jl) = P_v(jl) / P_id;

        %a2_v_0 = a2_v;

    end

    %a1_v_0 = a1_v;
    %a2_v_0 = a2_v;

    figure(1)
    plot(a1_v,'r')
    hold on
    plot(a2_v,'m')
    grid on
    title('a_D_1 & a_D_2 for each stream tube')
    xlabel('N of stream tube')
    ylabel('a_D')
    legend('a_D_1', 'a_D_2', 'location', 'Northeast')
    fprintf('a_D_1 and a_D_2 for V0 = %f and lam = %d \n', V0, lam)

    figure(2)
    plot(lam_v, Cp_v, '-o')
    hold on
    grid on
    ylim([0 1])
    xlim([lam_min lam_Max])
    xlabel('lambda [-]')
    ylabel('Cp [-]')
    legend('V0=5m/s','V0=7m/s','V0=9m/s','V0=11m/s','V0=13m/s','V0=15m/s','Location','northeast')
%%
    figure(3)
    plot(lam_v, P_v, '-o')
    hold on
    grid on
    %ylim([0 0.7])
    xlim([lam_min lam_Max])
    xlabel('lambda [-]')
    ylabel('P [kW]')
    legend('V0=5m/s','V0=7m/s','V0=9m/s','V0=11m/s','V0=13m/s','V0=15m/s','Location','northeast')

end

