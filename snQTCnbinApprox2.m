function y = snQTCnbinApprox2(s,Q,T,Kr,K0,L,lamda,pl,h,p)
y = snQTCnbin(s,Q,T,Kr,0,L,lamda,pl,h,p);

mtdem = -lamda.*T.*pl/((1-pl).*log(1-pl));  % itc20181209: the formula used to NOT have the T mul. factor

K0_cost_approx = min(mtdem./Q, 1).*K0./T;  % itc20181209: the formula used to NOT have the T div. factor

y = y + K0_cost_approx;

end
