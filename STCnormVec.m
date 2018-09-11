function y=STCnormVec(v,Kr,K0,L,mi,sigma,h,p,p2)
% return cost of (S,T) policy for normal demands
% v is the policy parameters vector [s T] 
s=v(1);
T=v(2);
Q=0;
P0=1; % Porder is always 1
Pi1=h.*(s+Q/2-mi.*(L+T/2));
Pi2=quadgk(@(t) staux(t,s,L,mi,sigma),0,T);
Pi2 = ((h+p)./(2*T)).*Pi2;
Pi3 = 0;
if nargin == 9 && p2 ~= 0
    Pi3 = quadgk(@(t) E(t,s,L+T,mi,sigma),s,Inf);
    Pi3 = p2*Pi3./T;
end
y = (Kr+K0*P0)./T + Pi1 + Pi2 + Pi3;
end


function z=staux(t,s,L,mi,sigma)
c = (s-mi.*(L+t))./(sigma.*sqrt(L+t));
der = 2.*(c.*(normcdf(c)-1)+normpdf(c));
z = sigma.*sqrt(L+t).*der;
end

function z=E(x,s,LT,mi,sigma)
z=(x-s).*normpdf((x-mi.*LT)./(sigma.*sqrt(LT)));
end

