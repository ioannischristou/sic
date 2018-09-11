function [y l u]=sSTCpoisson(s,S,T,Kr,K0,L,lamda,h,p,p2,eps)
% p2 is the p in the Hadley-Whittin notation
% p is the phat in the Hadley-Whittin notation
if nargin < 10
    p2 = 0;
end
% eps is the convergence accuracy
if nargin <11
    eps = 10^-8; % was 10^-9
end
% convert into Hadley-Whittin notation
J=Kr;
A=K0;
r=s;
R=S;
l=lamda;
m=lamda*L;
t=L;
phat=p;
IC=h;

K1 = J/T;
nom = sum1(r,R,T,t,l,m,IC,phat,p2,eps);
denom = sum2(r,R,T,l,eps);

y = K1 + (A+nom)/(T*denom);
l = K1 + nom/(T*denom);
u = K1 + A/T + nom/(T*denom);

end


% auxiliary functions


function y = sum1(r,R,T,t,l,m,IC,phat,p2,eps)
y = 0;
n=0;
last = 0; count = 0;
while true
    sum=0;
    for j=1:R-r
        sum = sum + poisspdf(R-r-j,n*l*T)*H(r+j,T,t,l,m,IC,phat,p2);
    end
    y = y+sum;
    last = last+sum;
    count = count + 1;
    if count == 4
        if y==0 && last==0
            break;
        end
        err = abs(last/y);
        if err < eps
            break;  % series converged
        else 
            count = 0;
            last = 0;
        end
    end
    n = n+1;
end
%disp(['sum1() took ' num2str(n) ' steps to converge.']);
end


function y = H(rpj,T,t,l,m,IC,phat,p2)
z = 0;
if nargin == 8 && p2~=0
    z = p2*eP(rpj,T,l,t);
end
y = IC*T*(rpj - m - l*T/2) + (IC+phat)*bP(rpj,T,l,t) + z;
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


function y = sum2(r,R,T,l,eps)
y = 0;
n=1;
last = 0; count = 0;
while true
    sum=0;
    for j=1:R-r
        sum = sum + n*poisspdf(R-r-j,(n-1)*l*T)*complpoisscdf(j,l*T);
    end
    y = y+sum;
    last = last+sum;
    count = count + 1;
    if count == 4
        if abs(last/y) < eps
            break;  % series converged
        else 
            count = 0;
            last = 0;
        end
    end
    n = n+1;
end
%disp(['sum2() took ' num2str(n) ' steps to converge.']);
end


function y=complpoisscdf(x,lm)
y = 1.0 - poisscdf(x-1,lm);
end