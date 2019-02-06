% comparison between exact (s,S,T) and approximation algorithms

fid = fopen('snQTCnormApproxSA.m10.csv', 'wt');

L=1;
lamda = 10;
mi = lamda;
sigma = sqrt(lamda);
%tnot = 10^-6;
tnot = (3.5*sigma/mi)^2;
p2 = 0;
Tmax = 5.0;

kmax = 10000;  % total number of function evaluations in SA

% Kr K0 h p s* S* T* c* time* 

fprintf(fid, '    Kr ,     Ko ,      h ,      p ,    s ,    S ,    T ,    c , time\n');

% Krs must be in descending order!
Krs=[3 2 0];
Kr_len = numel(Krs);
K0s=[5 25 100 1000];
K0_len = numel(K0s);
hs=[10 15 20];
h_len = numel(hs);
ps = [20 25 100];
p_len = numel(ps);

for j=1:K0_len
    K0=K0s(j);
    for k=1:h_len
        h=hs(k);
        for l=1:p_len
            p=ps(l);
            if p==100 && h~=15
                continue;  % only run h=15,p=100 combinations
            end
            if p==20 && h~=20
                continue;
            end
            if p==25 && h~=10
                continue;
            end
            x0=[0; L*lamda; 1];
            for i=1:Kr_len
                Kr=Krs(i);
                tic;
                [xb yb xs k] = nonlinoptSA(@(v) snQTCnorm_2(v,Kr,K0,L,mi,sigma,h,p,p2), x0, kmax, 1, @(x,h) quantrandgen(x,h));
                timesi=toc;
                si = xb(1,1);
                Qi = xb(2,1);
                Ti = xb(3,1);
                % bring into the feasible region
                if Ti<tnot
                    Ti=tnot;
                end
                if Qi<0
                    Qi=0;
                end
                ci=yb;
                fprintf(fid,'%6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f ,', Kr, K0, h, p, si, Qi, Ti, ci, timesi);
                fprintf(fid,'\n');
            end
        end
    end
end

fclose(fid);


function y=snQTCnorm_2(v,Kr,K0,L,mi,sigma,h,p,p2)
r=v(1);
Q=v(2);
if Q<0
    Q=1.e-8;
end
T=v(3);
tnot = (3.5*sigma/mi)^2;
if T<tnot
    T=tnot;
end
y = snQTCnorm(r,Q,T,Kr,K0,L,mi,sigma,h,p,p2);
end


function d = quantrandgen(x, h)
n = numel(x);
d = zeros(n,1);
epsQ = 0.1;
epsT = 0.01;
for i=1:n
    d(i,1) = random('unif',-h,h);
end
d(2,1) = round(d(2,1)/epsQ)*epsQ;
d(3,1) = round(d(3,1)/epsT)*epsT;
end
