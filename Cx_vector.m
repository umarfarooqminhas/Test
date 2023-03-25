function [a_v, Cx_v] = Cx_vector (Cx_a,at,step_a)

a_v_1 = 0:step_a:at;
a_v_2 = at+step_a:step_a:1;

for j1 = 1:length(a_v_1)
    a = a_v_1(j1);
    cx_1(j1) = 4*(a*(1-a));
end

for j2 = 1:length(a_v_2)
    a = a_v_2(j2);
    cx_2(j2) = Cx_a - 4*(1-a)*((Cx_a^(0.5))-1);
end

a_v = [a_v_1,a_v_2]';

Cx_v = [cx_1,cx_2]';


