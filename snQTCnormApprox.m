function y=snQTCnormApprox(s,Q,T,Kr,K0,L,mi,sigma,h,p)
% (nQ,r,T) model in continuous time T
% uses the Hadley-Whittin formula 5-63.
% and an Approximation formula for P0=Prob{Order}

% P0
% RT=Q/(sigma*sqrt(T));
% MT=T*mi/(sigma*sqrt(T));
% P0=1-(1/RT)*(normpdf(RT-MT)+(RT-MT)*normcdf(RT-MT)-(normpdf(-MT)-MT*normcdf(-MT)));
P0 = mi*T/Q;
if P0 > 1
    P0 = 1;
end

% Holding costs
H1 = h*(Q/2 + s - L*mi - mi*T/2);

% backorder costs
K1 = KsiA(s, L+T, mi, sigma^2);
K2 = KsiA(s, L, mi, sigma^2);
K3 = KsiA(s+Q, L+T, mi, sigma^2);
K4 = KsiA(s+Q, L, mi, sigma^2);
Yr = (K1 - K2)/T;
YrQ = (K3 - K4)/T;
B = (Yr - YrQ)/Q;
V1 = lamdaA(s,L+T, mi, sigma^2);
V2 = lamdaA(s,L, mi, sigma^2);
V3 = lamdaA(s+Q,L+T, mi, sigma^2);
V4 = lamdaA(s+Q,L, mi, sigma^2);
hr = (V1 - V2)/T;
hrQ = (V3 - V4)/T;
E = (hr - hrQ)/Q;
% B1 = h*B + p*E;
B1 = (h+p)*B;
%total costs
y = (Kr+K0*P0)/T + H1 + B1;
end

function z = KsiA(r, tau, l, D)
z1 = (1-normcdf((r-l*tau)/sqrt(D*tau)))*(l^2*tau^3/6 - D^2*r/(4*l^3) - l*tau^2*r/2 - D*r^2/(4*l^2) + D*tau^2/4 + tau*r^2/2 - r^3/(6*l) - D^3/(8*l^4));
z2 = normpdf((r-l*tau)/sqrt(D*tau))*(l*sqrt(D)*(sqrt(tau))^5/6 - sqrt(D)*(sqrt(tau))^3*r/3 + sqrt(D)*sqrt(tau)*r^2/(6*l) + (sqrt(D))^3*(sqrt(tau))^3/(12*l) + (sqrt(D))^3*sqrt(tau)*r/(4*l^2) + (sqrt(D))^5*sqrt(tau)/(4*l^3));
z3 = (1-normcdf((r+l*tau)/sqrt(D*tau)))*(D^3*exp(2*l*r/D)/(8*l^4));
z = z1+z2+z3;
end

function w = lamdaA(v, t, l, D)
waux = (1+((v-l*t)/sqrt(D*t))^2)*(1-normcdf((v-l*t)/sqrt(D*t))) - ((v-l*t)/sqrt(D*t))*normpdf((v-l*t)/sqrt(D*t));
w = waux*D*t/2;
end