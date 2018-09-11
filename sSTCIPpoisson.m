function y = sSTCIPpoisson(rpj,T,t,l,m,IC,phat,p2)
z = 0;
if nargin == 8 && p2~=0
    z = p2*eP(rpj,T,l,t);
end
y = IC*T*(rpj - m - l*T/2) + (IC+phat)*bP(rpj,T,l,t) + z;
y=y/T; % itc: HERE didn't use to have this line
end


function y = bP(rpj,T,l,t)
t1  = l*(((t+T)^2)*complpoisscdf(rpj-1,l*(t+T))-(t^2)*complpoisscdf(rpj-1,l*t))/2;
t2 = (complpoisscdf(rpj+1,l*(t+T))-complpoisscdf(rpj+1,l*t))*rpj*(rpj+1)/(2*l);
t3 = rpj*((t+T)*complpoisscdf(rpj,l*(t+T))-t*complpoisscdf(rpj,l*t));
y = t1 + t2 - t3;
end

function y = eP(rpj,T,l,t)
t1 = l*(t+T)*complpoisscdf(rpj-1,l*(t+T)) - rpj*complpoisscdf(rpj,l*(t+T));
t2 = l*t*complpoisscdf(rpj-1,l*t) - rpj*complpoisscdf(rpj,l*t);
y = t1 - t2;
end

function y=complpoisscdf(x,lm)
y = 1.0 - poisscdf(x-1,lm);
end