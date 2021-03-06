function y=snQTCnbin(s,Q,T,Kr,K0,L,lamda,p,h,phat)
% (nQ,r,T) model in continuous time T
% uses the Hadley-Whittin formulas 5-20, 5-22...24
% 5-27...29, and 5-33.
pm = 1- p;
r = -lamda.*L./log(1-p);
%P0 = lamda*T*nbincdf(Q-1,-lamda*T/log(1-p),pm)/Q + complnbincdf(Q+1,-lamda*T/log(1-p),pm);
P0=porder(Q,T,lamda,p);
ltdem = r.*(1-pm)./pm;  % mean lead-time demand
%y1 = (Kr+K0*P0)./T + h*((Q+1)/2.0  + s - ltdem - lamda*T/2.0);
% itc20181023: y1 used to be as above, but the last term doesn't capture 
% period mean-demand
perdem = -lamda.*T.*p/((1-p).*log(1-p));
y1 = (Kr+K0.*P0)./T + h.*((Q+1)./2.0  + s - ltdem - perdem./2.0);

% y = y1 + p*EP(s,Q,T,L,lamda) + h*BP(s,Q,T,L,lamda);
B = (h+phat).*BP(s,Q,T,L,lamda,p);
%disp(['y1=' num2str(y1) ' B=' num2str(B)]);
y = y1 + B;

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


function z1=BP(s,Q,T,L,lamda,p)
intgl = 0;
sp1 = s+1;
spQ = s+Q+1.e-9;  % itc20181209: added eps in right-end of sum
for u=sp1:spQ
    ipu = IP(T,L,lamda,p,u);
    %disp(['ipu(' num2str(u) ')=' num2str(ipu)]);  % itc: HERE rm asap
    intgl = intgl + ipu;
end
z1 = intgl./(Q.*T);
end


% function z = IP(T,L,lamda,p,u)
% z = quadgk(@(t) IIP(t,lamda,p,u), L, L+T, 'AbsTol', 1.e-12, 'RelTol', 0, 'MaxIntervalCount', 2000);  % itc20181209: integration params added
% %z = quadvec(@(t) IIP(t,lamda,p,u), L, L+T);
% end

function z = IP(T,L,lamda,p,u)
%z = quadgk(@(t) IIP(t,lamda,p,u), L, L+T, 'AbsTol', 1.e-12, 'RelTol', 0, 'MaxIntervalCount', 2000);  % itc20181209: integration params added
% itc20181210: break the integration interval into many sub-intervals...
num_intvls = 5;
len = T/num_intvls;
z=0;
left=L;
for i=1:num_intvls
    right = left+len;
    z = z + quadgk(@(t) IIP(t,lamda,p,u), left, right);
    left=right;
end
end


function z = IIP(t,lamda,p,u)
pm = 1-p;
r = -lamda.*t./log(1-p);
mean = p.*r./(1-p);
sum=0;
u1 = u-1+1.e-9;  % itc20181209: added eps in right-end of sum
for x=1:u1
    sum = sum + x.*nbinpdf(x,r,pm); 
    % itc20181002: last arg used to be pm, now is p
    % itc20181023: argument is correctly pm, because MATLAB computes for
    % nbin the function (r+k-1)C(k) p^r (1-p)^k thus for MATLAB nbinpdf/cdf
    % the p parameter is really 1-p.
end
ms = mean - sum;
last_term = u.*complnbincdf(u,r,p); % itc20181002: last arg used to be pm, now is p
z = ms - last_term; % itc20181002: last arg used to be pm, now is p
% if z<-1.e-6
%     disp(['IIP(t=' num2str(t) ',lambda=' num2str(lamda) ',p=' num2str(p) ',u=' num2str(u) ')=' num2str(z) ' (mean=' num2str(mean) ' sum=' num2str(sum) ' last_term=' num2str(last_term) ')']);  % itc: HERE rm asap
%     error('negative IIP');
% end
%disp(['IIP(' num2str(u) ')=' num2str(z) ' (mean=' num2str(mean) ' sum=' num2str(sum) ')']);  % itc: HERE rm asap
end

function y=complnbincdf(x,r,p)
pm = 1-p;
%if x==0
%    y=0;
%else
    y = 1.0 - nbincdf(x-1,r,pm);
%end
end
