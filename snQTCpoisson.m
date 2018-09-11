function y=snQTCpoisson(s,Q,T,Kr,K0,L,lamda,h,phat,p)
% (nQ,r,T) model in continuous time T
% uses the Hadley-Whittin formulas 5-20, 5-22...24
% 5-27...29, and 5-33.

if nargin < 10
    p = 0;
end

P0 = lamda*T*poisscdf(Q-1,lamda*T)/Q + complpoisscdf(Q+1,lamda*T);

y1 = (Kr+K0*P0)./T + h*((Q+1)/2.0  + s - lamda*L - lamda*T/2.0);
% y = y1 + p*EP(s,Q,T,L,lamda) + h*BP(s,Q,T,L,lamda);
y = y1 + (h+phat)*BP(s,Q,T,L,lamda) + p*EP(s,Q,T,L,lamda);

end


function z1=BP(s,Q,T,L,lamda)
z1 = (YP(s,T,L,lamda)-YP(s+Q,T,L,lamda))./Q;
end

function z2=YP(v,T,L,lamda)
z2 = KsiP(v,T+L,lamda,T)-KsiP(v,L,lamda,T);
end

function z3=KsiP(v,t,lamda,Tcap)
z31 = -lamda*v*(t^2)*complpoisscdf(v,lamda*t)./(2*Tcap);
z32 = (lamda^2)*(t^3)*complpoisscdf(v-1,lamda*t)./(6*Tcap);
z33 = v*(v+1)*t*complpoisscdf(v+1,lamda*t)./(2*Tcap);
z34 = - v*(v+1)*(v+2)*complpoisscdf(v+2,lamda*t)./(6*lamda*Tcap);
z3 = z31+z32+z33+z34;
end

function w=EP(s,Q,T,L,lamda)
w = (LamdaP(s,T,L,lamda) - LamdaP(s+Q,T,L,lamda))./Q;
end

function w2=LamdaP(v,T,L,lamda)
w2 = (betaP(v,L+T,lamda)-betaP(v,L,lamda))./T;
end

function w3=betaP(v,t,lamda)
w31 = ((lamda*t).^2).*complpoisscdf(v-1,lamda.*t)/2.0;
w32 = lamda.*t.*v.*complpoisscdf(v,lamda.*t);
w33 = v.*(v+1).*complpoisscdf(v+1,lamda.*t)./2.0;
w3 = w31 - w32 + w33;
end

function y=complpoisscdf(x,lm)
y = 1.0 - poisscdf(x-1,lm);
end

