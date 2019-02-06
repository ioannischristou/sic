function [c hm]=snQTCnbin2(s,Q,T,Kr,K0,L,lamda,p,h,phat)
c = snQTCnbin(s,Q,T,Kr,K0,L,lamda,p,h,phat);
P0=porder(Q,T,lamda,p);
hm = c - (Kr+K0*P0)./T;
end

function z = porder(Q,T,lamda,p)
%pm = 1-p;
z = 0;
Qp = Q+1.e-9;
for j=1:Qp
    z = z+complnbincdf(j,-lamda.*T./log(1-p),p);  % itc20181002: last arg used to be pm, now is p
end
z = z./Q;
end


function y=complnbincdf(x,r,p)
pm = 1-p;
if x==0
    y=0;
else
    y = 1.0 - nbincdf(x-1,r,pm);  
    % itc20181210: see snQTCnbin for an explanation why the last argument
    % must be pm instead of p
end
end


