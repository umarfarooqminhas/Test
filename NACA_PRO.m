function [av,lv,dv,Ev,Emax] = NACA_PRO (Re_i, NACA_40K, NACA_80K, NACA_160K, NACA_350K, NACA_700K, NACA_1M, NACA_2M, NACA_5M)

Re = Re_i;

av = [];
lv = [];
dv = [];
Ev = [];
Emax=0;

if Re < 40000    % Re lower than the lower reynolds available
    M1 = NACA_40K;
    Re1 = 40000;
    M2 = NACA_80K;
    Re2 = 80000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re );

elseif Re == 40000                % Re = 40k
    av = NACA_40K(:,1);
    lv = NACA_40K(:,2);
    dv = NACA_40K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);

elseif Re > 40000 && Re < 80000   % Re btwn 40K and 80K
    M1 = NACA_40K;
    Re1 = 40000;
    M2 = NACA_80K;
    Re2 = 80000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re );

elseif Re == 80000                % Re = 80k
    av = NACA_80K(:,1);
    lv = NACA_80K(:,2);
    dv = NACA_80K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);

elseif Re > 80000 && Re < 160000  % Re btwn 80k and 160k
    M1 = NACA_80K;
    Re1 = 80000;
    M2 = NACA_160K;
    Re2 = 160000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re ); 

elseif Re == 160000               % Re = 160k
    av = NACA_160K(:,1);
    lv = NACA_160K(:,2);
    dv = NACA_160K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);

elseif Re > 160000 && Re < 350000  % Re btwn 160k and 350k
    M1 = NACA_160K;
    Re1 = 160000;
    M2 = NACA_350K;
    Re2 = 350000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re ); 

elseif Re == 350000               % Re = 350k
    av = NACA_350K(:,1);
    lv = NACA_350K(:,2);
    dv = NACA_350K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);
  
elseif Re > 350000 && Re < 700000  % Re btwn 350k and 700K
    M1 = NACA_350K;
    Re1 = 350000;
    M2 = NACA_700K;
    Re2 = 700000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re ); 

elseif Re == 700000               % Re = 700k
    av = NACA_700K(:,1);
    lv = NACA_700K(:,2);
    dv = NACA_700K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);

elseif Re > 700000 && Re < 1000000  % Re btwn 700K and 1M  
    M1 = NACA_700K;
    Re1 = 700000;
    M2 = NACA_1M;
    Re2 = 1000000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re ); 
    
elseif Re == 1000000               % Re = 1M  
    av = NACA_1M(:,1);
    lv = NACA_1M(:,2);
    dv = NACA_1M(:,3);
    Ev = lv./dv;
    Emax = max(Ev);

elseif Re > 1000000 && Re < 2000000   % Re btwn 1M and 2M  
    M1 = NACA_1M;
    Re1 = 1000000;
    M2 = NACA_2M;
    Re2 = 2000000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re );

elseif Re == 2000000               % Re = 2M     
    av = NACA_2M(:,1);
    lv = NACA_2M(:,2);
    dv = NACA_2M(:,3);
    Ev = lv./dv;
    Emax = max(Ev);    

elseif Re > 2000000 && Re < 5000000   % Re btwn 2M and 5M   
    M1 = NACA_2M;
    Re1 = 2000000;
    M2 = NACA_5M;
    Re2 = 5000000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re );    

elseif Re >= 5000000    % Re equal or higher than the higher reynolds available
    av = NACA_5M(:,1);
    lv = NACA_5M(:,2);
    dv = NACA_5M(:,3);
    Ev = lv./dv;
    Emax = max(Ev);
end
