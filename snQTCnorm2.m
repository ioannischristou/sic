function [y hm] = snQTCnorm2(s,Q,T,Kr,K0,L,mi,sigma,h,p)

if K0==0  % short-cut
    P0=0;
else if Q>1.e-6
%     RT=Q./(sigma.*sqrt(T));
%     MT=T.*mi./(sigma.*sqrt(T));
%     P0=1-(1/RT).*(normpdf(RT-MT)+(RT-MT).*normcdf(RT-MT)-(normpdf(-MT)-MT.*normcdf(-MT)));
% itc20181214: computation of P0 according to Hadley-Wittin 1963:
      qltsrt = (Q-mi.*T)./(sigma.*sqrt(T));
      npdf = normpdf(qltsrt);
      ncdf = normcdf(qltsrt);
      P0 = mi.*T.*ncdf./Q + (1-ncdf) - sigma.*sqrt(T).*npdf./Q;
% itc20181214: end computatio nof P0
    else 
        P0=1;  % Q <= 1.e-6 is meant for (R,T) policy, Q=0
    end
end
% Pi1=h*(s+Q/2-mi*(L+(T+1)/2));
Pi1 = h.*(s+Q/2-mi*(L+T/2));

if Q>1.e-6
    Pi2 = quadgk(@(t) penint(t,s,Q,L,mi,sigma), 0, T);  % was quad
    Pi2 = ((h+p)./(2*T)).*Pi2;  
    y = (Kr+K0.*P0)./T + Pi1 + Pi2;
    hm = y - K0.*P0./T;
else  % Q<=1.e-6 implies use the (R,T) policy
    y = (Kr+K0)./T + h.*(s-L.*mi-mi.*T./2) + (h+p).*B(s,T,L,mi,sigma);
    hm = y - K0./T;
end
end


function z = penint(T,s,Q,L,mi,sigma)
ZLpt = (s-mi*(L+T))./(sigma*sqrt(L+T));
RLpt = Q./(sigma*sqrt(L+T));
t1=((ZLpt+RLpt).^2+1).*normcdf(ZLpt+RLpt)+(ZLpt+RLpt).*normpdf(ZLpt+RLpt);
t2=(ZLpt.^2+1).*normcdf(ZLpt)+ZLpt.*normpdf(ZLpt)+(2*ZLpt+RLpt).*RLpt;
z = (t1 - t2).*sigma.*sqrt(L+T)./RLpt; 
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
% values of Q in the function snQTCnorm2(), so I prefer to keep it as is
end

