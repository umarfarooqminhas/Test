function [av,lv,dv,Ev,Emax] = DU_PRO (Re_i, DU_40K, DU_80K, DU_160K, DU_360K, DU_700K, DU_1M, DU_2M, DU_5M)

Re = Re_i;

av = [];
lv = [];
dv = [];
Ev = [];
Emax=0;

if Re < 40000    % Re lower than the lower reynolds available
    M1 = DU_40K;
    Re1 = 40000;
    M2 = DU_80K;
    Re2 = 80000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re );

elseif Re == 40000                % Re = 40k
    av = DU_40K(:,1);
    lv = DU_40K(:,2);
    dv = DU_40K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);

elseif Re > 40000 && Re < 80000   % Re btwn 40K and 80K
    M1 = DU_40K;
    Re1 = 40000;
    M2 = DU_80K;
    Re2 = 80000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re );

elseif Re == 80000                % Re = 80k
    av = DU_80K(:,1);
    lv = DU_80K(:,2);
    dv = DU_80K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);

elseif Re > 80000 && Re < 160000  % Re btwn 80k and 160k
    M1 = DU_80K;
    Re1 = 80000;
    M2 = DU_160K;
    Re2 = 160000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re ); 

elseif Re == 160000               % Re = 160k
    av = DU_160K(:,1);
    lv = DU_160K(:,2);
    dv = DU_160K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);

elseif Re > 160000 && Re < 360000  % Re btwn 160k and 360k
    M1 = DU_160K;
    Re1 = 160000;
    M2 = DU_360K;
    Re2 = 360000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re ); 

elseif Re == 360000               % Re = 360k
    av = DU_360K(:,1);
    lv = DU_360K(:,2);
    dv = DU_360K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);
  
elseif Re > 360000 && Re < 700000  % Re btwn 360k and 700K
    M1 = DU_360K;
    Re1 = 360000;
    M2 = DU_700K;
    Re2 = 700000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re ); 

elseif Re == 700000               % Re = 700k
    av = DU_700K(:,1);
    lv = DU_700K(:,2);
    dv = DU_700K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);

elseif Re > 700000 && Re < 1000000  % Re btwn 700K and 1M  
    M1 = DU_700K;
    Re1 = 700000;
    M2 = DU_1M;
    Re2 = 1000000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re ); 
    
elseif Re == 1000000               % Re = 1M  
    av = DU_1M(:,1);
    lv = DU_1M(:,2);
    dv = DU_1M(:,3);
    Ev = lv./dv;
    Emax = max(Ev);

elseif Re > 1000000 && Re < 2000000   % Re btwn 1M and 2M  
    M1 = DU_1M;
    Re1 = 1000000;
    M2 = DU_2M;
    Re2 = 2000000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re );

elseif Re == 2000000               % Re = 2M     
    av = DU_2M(:,1);
    lv = DU_2M(:,2);
    dv = DU_2M(:,3);
    Ev = lv./dv;
    Emax = max(Ev);    

elseif Re > 2000000 && Re < 5000000   % Re btwn 2M and 5M   
    M1 = DU_2M;
    Re1 = 2000000;
    M2 = DU_5M;
    Re2 = 5000000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re );    

elseif Re >= 5000000    % Re equal or higher than the higher reynolds available
    av = DU_5M(:,1);
    lv = DU_5M(:,2);
    dv = DU_5M(:,3);
    Ev = lv./dv;
    Emax = max(Ev);
end
