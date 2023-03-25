function [P_Re_1,P_Re_2,P_Re_3,P_Re_4,P_Re_5,P_Re_6,P_Re_7,P_Re_8] = extractionV (P_Re_1,P_Re_2,P_Re_3,P_Re_4,P_Re_5,P_Re_6,P_Re_7,P_Re_8)

s = 0.05;
a_first = -180;
a_end = 180;

[P_Re_1] = palas(P_Re_1,s, a_first, a_end);
[P_Re_2] = palas(P_Re_2,s, a_first, a_end);
[P_Re_3] = palas(P_Re_3,s, a_first, a_end);
[P_Re_4] = palas(P_Re_4,s, a_first, a_end);
[P_Re_5] = palas(P_Re_5,s, a_first, a_end);
[P_Re_6] = palas(P_Re_6,s, a_first, a_end);
[P_Re_7] = palas(P_Re_7,s, a_first, a_end);
[P_Re_8] = palas(P_Re_8,s, a_first, a_end);