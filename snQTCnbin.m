function y=snQTCnbin(s,Q,T,Kr,K0,L,lamda,p,h,phat)
% (nQ,r,T) model in continuous time T
% uses the Hadley-Whittin formulas 5-20, 5-22...24
% 5-27...29, and 5-33.
pm = 1- p;
r = -lamda*L/log(1-p);
%P0 = lamda*T*nbincdf(Q-1,-lamda*T/log(1-p),pm)/Q + complnbincdf(Q+1,-lamda*T/log(1-p),pm);
P0=porder(Q,T,lamda,p);
ltdem = r*(1-pm)/pm;  % mean lead-time demand
y1 = (Kr+K0*P0)./T + h*((Q+1)/2.0  + s - ltdem - lamda*T/2.0);
% y = y1 + p*EP(s,Q,T,L,lamda) + h*BP(s,Q,T,L,lamda);
B = (h+phat)*BP(s,Q,T,L,lamda,p);
%disp(['y1=' num2str(y1) ' B=' num2str(B)]);
y = y1 + B;

end


function z = porder(Q,T,lamda,p)
pm = 1-p;
z = 0;
Qp = Q+1.e-9;
for j=1:Qp
    z = z+complnbincdf(j,-lamda*T/log(1-p),pm);
end
z = z/Q;
end


function z1=BP(s,Q,T,L,lamda,p)
intgl = 0;
for u=s+1:s+Q
    ipu = IP(T,L,lamda,p,u);
    %disp(['ipu(' num2str(u) ')=' num2str(ipu)]);  % itc: HERE rm asap
    intgl = intgl + ipu;
end
z1 = intgl/(Q*T);
end


function z = IP(T,L,lamda,p,u)
z = quadgk(@(t) IIP(t,lamda,p,u), L, L+T);
end

function z = IIP(t,lamda,p,u)
pm = 1-p;
r = -lamda*t/log(1-p);
mean = p*r/(1-p);
sum=0;
u1 = u-1;
for x=1:u1
    sum = sum + x*nbinpdf(x,r,pm);
end
ms = mean - sum;
z = ms - u*complnbincdf(u,r,pm);
%disp(['IIP(' num2str(u) ')=' num2str(z)]);  % itc: HERE rm asap
end

function y=complnbincdf(x,r,p)
y = 1.0 - nbincdf(x-1,r,p);
end

