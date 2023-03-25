function [a1_i_v_up,T_D1,TAO_1] = VAWT_UP_NACA (N_st, D_theta, theta_v_rad_up, a_v, Cx_v, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M)

tol_a = 0.001;
tol_alp = 0.04; %tol_cx = 0.05;  %0.005 %0.5 => find u

a1_i_old = 0.2;
a1_i = a1_i_old;
a1_i_v_up = a1_i;

e = 0.5; %0.8; %relaxation

it_max = 1000;

for i = 1 : N_st

    theta_i = theta_v_rad_up(i);
    theta= theta_i;

    dA_i = H * R * abs(sin(theta_i))* D_theta;

    it = 0;
    err_a = tol_a + 1;     

    %a1_i = a1_i_v(i);

    while it < it_max && err_a > tol_a

        V_D1_i = V0 * (1-a1_i);

        w_inf_t = V_D1_i*cos(theta_i)+omega*R; 
        w_inf_n = V_D1_i*sin(theta_i);
        w_D1_i = (w_inf_t^2+w_inf_n^2)^(1/2);
        w = w_D1_i;

        Re_i = w_D1_i *c*ro/mi;

        fi_D1_i = atan(((1-a1_i)*sin(theta_i))/((1-a1_i)*cos(theta_i)+lam));
        alfa_D1_i = fi_D1_i - 0;   %bc = 0 % radians
        alfa = rad2deg(alfa_D1_i);

        [av,lv,dv,Ev,Emax] = NACA_PRO (Re_i, NACA_40K, NACA_80K, NACA_160K, NACA_350K, NACA_700K, NACA_1M, NACA_2M, NACA_5M);
        %[av,lv,dv,Ev,Emax] = DU_PRO (Re_i, P_Re_1,P_Re_2,P_Re_3,P_Re_4,P_Re_5,P_Re_6,P_Re_7,P_Re_8);
        
        [k] = find_tol (av,alfa,tol_alp); %alpha is in degrees
        cl = lv(k);
        cd = dv(k);

        % finite blade correction
        alfa = alfa - rad2deg(abs(cl)/(pi*AR));
        %alfa = alfa - abs(cl)/(pi*AR);
                
        [kp] = find_tol (av,alfa,tol_alp);
        cl = lv(kp);
        cd = cd + (cl^2)/(pi*AR);

        Cn_D1_i = cl * cosd(alfa) + cd * sind(alfa);
        Ct_D1_i = cl * sind(alfa) - cd * cosd(alfa);

        dFn_D1_i = 0.5*ro*w^2*c*H*Cn_D1_i;
        dFt_D1_i = 0.5*ro*w^2*c*H*Ct_D1_i;

        dFx_D1_i = dFn_D1_i*sin(theta)-dFt_D1_i*cos(theta);
        dFx_D1_i_avg = (Nb/(2*pi))*D_theta*dFx_D1_i;

        Cx_D1_i = dFx_D1_i_avg / (0.5*ro*V0^2*dA_i);

        %[k_x] = find_tol(Cx_v,Cx_D1_i,tol_cx);
        [k_x] = find_U(Cx_v,Cx_D1_i);

        a1_i_new = a_v(k_x);

        
        if isnan(a1_i_new) 
            a1_i_new = a1_i_v_up(i-2);
        elseif a1_i_new > 1
            a1_i_new = a1_i_v_up(i-2);
        elseif a1_i_new < -1
            a1_i_new = a1_i_v_up(i-2);
        end

        %a1_i_new

        err_a = abs(a1_i_new - a1_i);
        
        it = it+1;
        a1_i = (e*a1_i_new+(1-e)*a1_i);
        %a1_i = a1_i_new;

    end

    %a1_i_v(i) = a1_i_new;
    a1_i_v_up(i) = a1_i_new;

    T_D1(i) = 0.5 * ro * w^2 *c*H*R* Ct_D1_i;
    A_D1(i) = dA_i;
    TAO_1(i) = T_D1(i) * A_D1(i) * omega;

end
