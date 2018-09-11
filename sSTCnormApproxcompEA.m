% comparison between exact (s,S,T) and approximation algorithms

fid = fopen('sSTCnormApproxEA.txt', 'wt');

L=1;
m = 50;
s = sqrt(50);
%tnot = 10^-6;
tnot = (3*s/m)^2;
epst=0.01;
epsq=10^-3;
p2 = 0;
estop = 1;
Tmax = 5.0;

% Kr K0 h p s* S* T* c* time* 

fprintf(fid, '    Kr ,     Ko ,      h ,      p ,    s ,    S ,    T ,    c , time\n');

% Krs must be in descending order!
Krs=[5 2 1];
Kr_len = numel(Krs);
K0s=[1 5 25 100 1000];
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
            x0=[0; L*m; 1];
            sigma_ea = [1; 1; 1];
            for i=1:Kr_len
                Kr=Krs(i);
                tic;
                [xb yb xs k] = nonlinoptEA(@(v) sSTCnorm_2(v,Kr,K0,L,m,s,h,p,p2), x0, sigma_ea, 200);  % default k is 1000
                timesi=toc;
                si = xb(1,1);
                Si = xb(2,1);
                Ti = xb(3,1);
                % bring into the feasible region
                if Ti<tnot
                    Ti=tnot;
                end
                if Si-si<=0
                    si = Si-epsq;
                end
                ci=yb;
                fprintf(fid,'%6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f ,', Kr, K0, h, p, si, Si, Ti, ci, timesi);
                fprintf(fid,'\n');
            end
        end
    end
end

fclose(fid);
