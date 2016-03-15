function RorQ = ORc(R, Q)
% OR for RFC conceptors R, Q 

notR = 1 - R; notQ = 1 - Q;

RorQ = 1 - (ANDc(notR, notQ));

