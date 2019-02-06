L=1;
lamda=0.5;
p=0.8;

Tmin=11;
Tmax=12;
epst = 0.1;
Tlen = (Tmax-Tmin)/epst;

s=6;
Q=1;
T=11.7;
Ti=1:Tlen;
BPi=1:Tlen;

% ns=50;
% Xi=1:ns;
% IIPi=1:ns;
% eps=T/ns;
% for i=1:ns
%     Xi(i)=L+eps*(i-1);
%     IIPi(i)=IIP(i,lamda,p,s+Q);
% end
% hold on
% plot(Xi,IIPi);
% hold off

% Xi=1:20;
% IPi=1:20;
% for i=1:20
%     Xi(i)=i;
%     IPi(i)=IP(T,L,lamda,p,i);
% end
% plot(Xi,IPi);

% for i=1:Tlen
%     T = Tmin+(i-1)*epst;
%     Ti(i)=T;
%     BPi(i)=BP(s,Q,T,L,lamda,p);
% end
% hold on
% plot(Ti,BPi);
% hold off

% Qi=1:50;
% mqi=1:50;
% P0i=1:50;
% diffi=1:50;
% for i=1:50
%     P0i(i)=porder(i,T,lamda,p);
%     mqi(i)=min(-lamda*T*p/((1-p)*log(1-p)*i),1);
%     diffi(i)=mqi(i)-P0i(i);
% end
% hold on
% %plot(Qi,P0i,'k');
% plot(Qi,diffi);
% hold off
% hold on
% plot(Qi,mqi,'r');
% hold off

Xi=1:10;
Bi=1:10;
for i=1:10
    ni = (i-1)*100+1;
    Bi(i)=BP(s,Q,T,L,lamda,p,ni);
end
plot(Xi,Bi);


function z1=BP(s,Q,T,L,lamda,p,num_intvls)
intgl = 0;
sp1 = s+1;
spQ = s+Q+1.e-9;  % itc20181209: added eps in right-end of sum
for u=sp1:spQ
    ipu = IP(T,L,lamda,p,u,num_intvls);
    %disp(['ipu(' num2str(u) ')=' num2str(ipu)]);  % itc: HERE rm asap
    intgl = intgl + ipu;
end
z1 = intgl./(Q.*T);
end


function z = IP(T,L,lamda,p,u, num_intvls)
%z = quadgk(@(t) IIP(t,lamda,p,u), L, L+T, 'AbsTol', 1.e-12, 'RelTol', 0, 'MaxIntervalCount', 2000);  % itc20181209: integration params added
% itc20181210: break the integration interval into many sub-intervals...
if nargin < 6
    num_intvls = 5;
end
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


function z = porder(Q,T,lamda,p)
%pm = 1-p;
z = 0;
Qp = Q+1.e-9;
for j=1:Qp
    z = z+complnbincdf(j,-lamda.*T./log(1-p),p);  % itc20181002: last arg used to be pm, now is p
end
z = z./Q;
end

