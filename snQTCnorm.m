function [y, hm]=snQTCnorm(s,Q,T,Kr,K0,L,mi,sigma,h,p,p2,qnot)
% (nQ,r,T) model in continuous time T
% uses the Hadley-Whittin formula 5-63.
% Notice that when Q<=0.1 we assume Q=0, and we revert to computing the
% (R,T) formula according to Hadley-Whittin formula 5-74. In such a case,
% we assume that p2=0, for we don't implement the formula 5-72. Also notice
% that there is a grey area in computing (r,nQ,T) for small values of Q,
% since for small enough Q, the "negligible term"
% normpdf((Q(=0)-mi*T)/(sigma*sqrt(T)) or the term
% (1-normcdf((Q(=0)-mi*T)/(sigma*sqrt(T)) will not be sufficiently close to
% zero to offset the number 1/Q, which will skyrocket!
% y is the total cost
% hm is the cost when K0=0 (cost minus ordering costs)
if nargin < 12
    qnot = 1.e-6;  % what value of Q is considered zero --> (R,T)
end
if nargin < 11
    p2 = 0;
end
% P0
if Q<=qnot
    P0=1;
else
    RT=Q/(sigma*sqrt(T));
    MT=T*mi/(sigma*sqrt(T));
    P0=1-(1/RT)*(normpdf(RT-MT)+(RT-MT)*normcdf(RT-MT)-(normpdf(-MT)-MT*normcdf(-MT)));
    if P0>1
        error(['Q=' num2str(Q) ' T=' num2str(T) ' P0=' num2str(P0) '?']);
    end
end
% itc20181214: computation of P0 according to Hadley-Wittin 1963:
% if Q<=1.e-6
%     P0=1;
% else
%       qltsrt = (Q-mi.*T)./(sigma.*sqrt(T));
%       npdf = normpdf(qltsrt);
%       ncdf = normcdf(qltsrt);
%       P0 = mi.*T.*ncdf./Q + (1-ncdf) - sigma.*sqrt(T).*npdf./Q;
% end
%disp(['P0=' num2str(P0)]);  % itc20181214: HERE rm asap
% itc20181214: end computation nof P0

% Holding costs
H1 = h*(Q/2 + s - L*mi - mi*T/2);

if Q<=qnot
    % compute (R,T)
    if p2~=0
        error('p2 must be zero');
    end
% Hadley-Whitin formula below
%     UrLpT = U(s,L+T,mi,sigma^2);
%     UrT = U(s,L,mi,sigma^2);
%     B = (UrLpT-UrT)./T;
% (R,T) according to itc calculations
    B2 = B(s,T,L,mi,sigma);
    hm = Kr./T + H1 + (h+p).*B2;
    if isnan(hm)
        error('hm==NaN');
    end
    y = K0./T + hm;
    return;
end

% backorder costs
K1 = Ksi(s, L+T, mi, sigma^2);
if isnan(K1)
    error(['s=' num2str(s) ' q=' num2str(q) ' T=' num2str(T) ' L=' num2str(L) ' h=' num2str(h) ' p=' num2str(p) ' K1 is NaN']);
end
K2 = Ksi(s, L, mi, sigma^2);
if isnan(K2)
    error(['s=' num2str(s) ' q=' num2str(q) ' T=' num2str(T) ' L=' num2str(L) ' h=' num2str(h) ' p=' num2str(p) ' K2 is NaN']);
end
K3 = Ksi(s+Q, L+T, mi, sigma^2);
if isnan(K3)
    error(['s=' num2str(s) ' q=' num2str(q) ' T=' num2str(T) ' L=' num2str(L) ' h=' num2str(h) ' p=' num2str(p) ' K3 is NaN']);
end
K4 = Ksi(s+Q, L, mi, sigma^2);
if isnan(K4)
    error(['s=' num2str(s) ' q=' num2str(q) ' T=' num2str(T) ' L=' num2str(L) ' h=' num2str(h) ' p=' num2str(p) ' K4 is NaN']);
end
Yr = (K1 - K2)/T;
if isnan(Yr)
    error(['s=' num2str(s) ' q=' num2str(q) ' T=' num2str(T) ' L=' num2str(L) ' h=' num2str(h) ' p=' num2str(p) ' Yr is NaN']);
end
YrQ = (K3 - K4)/T;
if isnan(YrQ)
    error(['s=' num2str(s) ' q=' num2str(q) ' T=' num2str(T) ' L=' num2str(L) ' h=' num2str(h) ' p=' num2str(p) ' YrQ is NaN']);
end
B = (Yr - YrQ)/Q;
if isnan(B)
    error(['s=' num2str(s) ' q=' num2str(q) ' T=' num2str(T) ' L=' num2str(L) ' h=' num2str(h) ' p=' num2str(p) ' B is NaN']);
end
if p2>0
    V1 = lamda(s,L+T, mi, sigma^2);
    V2 = lamda(s,L, mi, sigma^2);
    V3 = lamda(s+Q,L+T, mi, sigma^2);
    V4 = lamda(s+Q,L, mi, sigma^2);
    hr = (V1 - V2)/T;
    hrQ = (V3 - V4)/T;
    E = (hr - hrQ)/Q;
    B2 = p2*E;
else
    B2=0;
end
% B1 = h*B + p*E;
B1 = (h+p)*B;

%total costs
hm = Kr/T + H1 + B1 + B2;
y = (Kr+K0*P0)/T + H1 + B1 + B2;
if isnan(y)
    error(['s=' num2str(s) ' q=' num2str(q) ' T=' num2str(T) ' L=' num2str(L) ' h=' num2str(h) ' p=' num2str(p) ' y is NaN']);
end
end

function z = Ksi(r, tau, l, D)
z1 = (1-normcdf((r-l*tau)/sqrt(D*tau)))*(l^2*tau^3/6 - D^2*r/(4*l^3) - l*tau^2*r/2 - D*r^2/(4*l^2) + D*tau^2/4 + tau*r^2/2 - r^3/(6*l) - D^3/(8*l^4));
z2 = normpdf((r-l*tau)/sqrt(D*tau))*(l*sqrt(D)*(sqrt(tau))^5/6 - sqrt(D)*(sqrt(tau))^3*r/3 + sqrt(D)*sqrt(tau)*r^2/(6*l) + (sqrt(D))^3*(sqrt(tau))^3/(12*l) + (sqrt(D))^3*sqrt(tau)*r/(4*l^2) + (sqrt(D))^5*sqrt(tau)/(4*l^3));
z3 = (1-normcdf((r+l*tau)/sqrt(D*tau)))*(D^3*exp(2*l*r/D)/(8*l^4));
z = z1+z2+z3;
end

function w = lamda(v, t, l, D)
waux = (1+((v-l*t)/sqrt(D*t))^2)*(1-normcdf((v-l*t)/sqrt(D*t))) - ((v-l*t)/sqrt(D*t))*normpdf((v-l*t)/sqrt(D*t));
w = waux*D*t/2;
end

function y = U(u,tau,l,D)
ultsrL = (u-l.*tau)./sqrt(D*tau);
t1 = (1-normcdf(ultsrL)).*((D^2+2*l^4*tau^2)/(4*l^3) + (D-2*l^2*tau)*u/(2*l^2) + u^2/(2*l));
t2 = 0.5*(sqrt(D)*tau^(1.5) - D^(1.5)*sqrt(tau)/l^2 - sqrt(D*tau)*u/l)*normpdf(ultsrL);
t3 = D^2*exp(2*l*u/D)*(1-normcdf((u+l*tau)/sqrt(D*tau)))/(4*l^3);
y = t1 + t2 - t3;
end

function y=B(r,T,L,m,s)
y = (Ksi_prime(r,L,m,s)-Ksi_prime(r,L+T,m,s))./T;
end

function z=Ksi_prime(r,L,m,s)
D=s.^2;
tau=L;
l=m;
P1 = l.^2*(tau.^3)./6 - D.^2*r./(4*l.^3) - l.*tau.^2*r./2 - D.*r.^2/(4*l.^2) + D.*(tau.^2)/4 + tau.*(r.^2)./2 - (r.^3)./(6*l) - (D.^3)./(8*l.^4);
P2 = l.*sqrt(D).*((sqrt(tau)).^5)./6 - sqrt(D).*((sqrt(tau)).^3).*r./3 + sqrt(D).*sqrt(tau).*(r.^2)./(6*l) + (sqrt(D)).^3*((sqrt(tau))^3)./(12*l) + (sqrt(D))^3*sqrt(tau).*r./(4*l.^2) + (sqrt(D)).^5*sqrt(tau)./(4*l.^3);

srL=s.*sqrt(L);
firmlm = normpdf((r-L.*m)./srL);
z = (1-normcdf((r-L*m)/srL))*(-s^4/(4*m^3) - m*L^2/2 - r*s^2/(2*m^2) + r*L - r^2/(2*m));
z = z + P1*(-firmlm/srL);
z = z + firmlm.*(-s*sqrt(L)^3/3 + srL*r/(3*m) + s^2*srL/(4*m^2));
z = z + P2*(firmlm*(L*m-r)/(L*D));
exp2mrs2 = exp(2*m*r/D);
z = z + (D^3./(8*m^4))*exp2mrs2*((1-normcdf((r+L*m)/srL))*2*m/D - normpdf((r+L*m)/srL)/srL);  
% itc20181210: first D-term, D above, in the Hadley-Wittin formula must be
% D^3, but if set to this value, the formula does not agree with small
% values of Q in the function snQTCnorm2()...
end

