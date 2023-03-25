function [a2_v_0,T_D2,TAO_2] = VAWT_DOWN_NACA_pp_0 (N_st, D_theta, theta_v_rad_down,a1_v_0, at, Cx_a, Nb, c,ro,mi,H,R,AR, V0, lam, omega, NACA_40K,NACA_80K,NACA_160K,NACA_350K,NACA_700K,NACA_1M,NACA_2M,NACA_5M)

tol_a = 0.001; %tol_alp = 0.04;  %tol_cx = 0.05;

a2_i_old = 0.08; % 0.08;
a2_i = a2_i_old;
a2_v_0 = a2_i;

e = 0.1; %relaxation

it_max = 1000; %N_st_tot = 3

for i = 1 : N_st

    %i

    theta_i = theta_v_rad_down(i);
    theta= theta_i;

    dA_i = H * R * abs(sin(theta_i))* D_theta;

    a1_i = a1_v_0(i);

    V_E_i = V0 * (1-2*a1_i);
    
    it = 0;
    err_a = tol_a + 1;     

    %a1_i = a1_i_v(i);

    while it < it_max && err_a > tol_a

        V_D2_i = V_E_i * (1-a2_i);

        w_inf_t = V_D2_i*cos(theta_i)+omega*R;
        w_inf_n = V_D2_i*sin(theta_i);
        w_D2_i = (w_inf_t^2+w_inf_n^2)^(1/2);
        w = w_D2_i;

        Re_i = w_D2_i *c*ro/mi;

        fi_D2_i = atan( ((1-a2_i)*sin(theta_i)) / ((1-a2_i)*cos(theta_i)+(lam/(1-2*a1_i))) );
        %fi_D2_i = atan( ((1-a2_i)*sin(theta_i)) / ((1-a2_i)*cos(theta_i)+(lam)) );

        alfa_D2_i = fi_D2_i - 0;   %bc = 0 % radians
        alfa = rad2deg(alfa_D2_i);

        [av,lv,dv,Ev,Emax] = NACA_PRO (Re_i, NACA_40K, NACA_80K, NACA_160K, NACA_350K, NACA_700K, NACA_1M, NACA_2M, NACA_5M);
        %[av,lv,dv,Ev,Emax] = DU_PRO (Re_i, DU_40K, DU_80K, DU_160K, DU_360K, DU_700K, DU_1M, DU_2M, DU_5M);
        
        %[k] = find_tol (av,alfa,tol_alp);
        [k] = find_U(av,alfa); %alpha is in degrees
        cl = lv(k);
        cd = dv(k);

        % finite blade correction
        alfa = alfa - rad2deg(abs(cl)/(pi*AR));         %alfa = alfa - abs(cl)/(pi*AR);
        
        %[kp] = find_tol (av,alfa,tol_alp);
        [kp] = find_U (av,alfa);
        cl = lv(kp);
        cd = cd + (cl^2)/(pi*AR);

        Cn_D2_i = cl * cosd(alfa) + cd * sind(alfa);
        Ct_D2_i = cl * sind(alfa) - cd * cosd(alfa);

        dFn_D2_i = 0.5*ro*w^2*c*H*Cn_D2_i;
        dFt_D2_i = 0.5*ro*w^2*c*H*Ct_D2_i;

        dFx_D2_i = dFn_D2_i*sin(theta)-dFt_D2_i*cos(theta);
        dFx_D2_i_avg = (Nb/(2*pi))*D_theta*dFx_D2_i;

        %Cx_D2_i = dFx_D2_i_avg / (0.5*ro*V_E_i^2*dA_i);
        %[Cx_v_down, D_Cx] = Cx_corr (a1_i,a_v_1,a_v_2,ac);
        %tol_cx = D_Cx*1000
        %[k_x] = find_tol(Cx_v,Cx_D2_i,tol_cx);

        %Cx_D2_i = dFx_D2_i_avg / (0.5*ro*V_E_i^2*dA_i);

        Cx_D2_i = dFx_D2_i_avg / (0.5*ro*V0^2*dA_i);

        if a2_i <= at
            %p(1)=(2+sqrt(4-4*Cx_D2_i))/4;
            p(1)=(2+sqrt(4-4*Cx_D2_i/(1-2*a1_i)^2))/4;
        
            %p(2)=(2-sqrt(4-4*Cx_D2_i))/4;
            p(2)=(2-sqrt(4-4*Cx_D2_i/(1-2*a1_i)^2))/4;

            p=p(p>0);
            a2_i_new = min(p);
        else
            %a2_i_new = 1-(Cx_a-Cx_D2_i)/(4*(sqrt(Cx_a)-1));
            a2_i_new = 1-(Cx_a-Cx_D2_i/(1-2*a1_i)^2)/(4*(sqrt(Cx_a)-1));
            
        end

         
        if isnan(a2_i_new) 
            a2_i_new = a2_v_0(i-2);
        elseif a2_i_new > 1
            a2_i_new = a2_v_0(i-2);
        elseif a2_i_new < -1
            a2_i_new = a2_v_0(i-2);
        end

        %a2_i_new

        err_a = abs(a2_i_new - a2_i);
        
        it = it+1;
        a2_i = (e*a2_i_new+(1-e)*a2_i);
        a2_i = real(a2_i);
        %a2_i = a2_i_new;
    end

    %a1_i_v(i) = a1_i_new;
    a2_v_0(i) = a2_i_new;

    %a2_i = a2_i_old;
    %if i <= 1
     %   a2_i = a2_i_old;
    %elseif i <= 5 && i >= 2
     %   a2_i = a2_i_v_down(i);
    %else
    %    a2_i = a2_i_v_down(i);
    %end

    a2_i = a2_i_new;

    a2_v_0(i) = a2_i_new;
    a2_v_down(i) = real(a2_v_0(i));
    a2_v_0(i) = a2_v_down(i);
    
    T_D2(i) = 0.5 * ro * w^2 *c*H*R* Ct_D2_i;
    A_D2(i) = dA_i;
    TAO_2(i) = T_D2(i) * A_D2(i) * omega;

end