% comparison between exact (s,S,T) and approximation algorithms

fid = fopen('sSTCnormExactHeurComp.new.txt', 'wt');

L=1;
mi = 50;
sigma = sqrt(50);
%tnot = 10^-6;
tnot = (3*sigma/mi)^2;
epsd = 0.5;
epst=0.01;
p2 = 0;
Tmax = 5.0;

epsq = 0.5;
qnot = 1.e-3;

% Kr K0 h p s* S* T* c* time* s1 S1 T1 c1 time1 r2 r2+Q2 T2 c2

fprintf(fid, '    Kr ,     Ko ,      h ,      p ,    s* ,    S* ,    T* ,    c* , time*, s1 , S1 , T1 , c1 , time1 , r2, r2+Q2, T2 , c2 , time2 \n');

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
            % first compute the (s,S,T) policies for all Krs
            [si Si Ti ci timesi] = sSTCnormOpt_MultiKr(Tmax, Krs, K0, L, mi, sigma, h, p, p2, epsd, epst, tnot);
            % now compute for each Kr the (s,nQ,T) and (R,T) policies
            for i=1:Kr_len
                Kr=Krs(i);
                tic;
                disp(['computing Approx policies for Kr=' num2str(Kr)]);
                [sopt Sopt Topt copt] = sSTCnormOptFastApprox2(Kr, K0, L, mi, sigma, h, p, p2);
                tElapsed=toc;
                if copt < ci(i)
                    % round-off error
                    si(i)=sopt;
                    Si(i)=Sopt;
                    Ti(i)=Topt;
                    ci(i)=copt;
                end
                fprintf(fid,'%6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f ,', Kr, K0, h, p, si(i), Si(i), Ti(i), ci(i), timesi(i));
                fprintf(fid,'%6.2f , %6.2f , %6.2f , %6.2f , %6.2f ,', sopt, Sopt, Topt, copt, tElapsed);
                disp(['computing opt. (s,nQ,T) policy']);
                tic;
                [r2 Q2 T2 c2] = snQTCnormOptFastApprox2(Kr,K0,L,mi,sigma,h,p);
                s2 = r2;
                if Q2 == 0
                    Q2 = epsq;
                end
                S2 = r2+Q2;
                c2 = sSTCnorm(s2,S2,T2,Kr,K0,L,mi,sigma,h,p);
                tElapsed = toc;
                fprintf(fid,'%6.2f , %6.2f , %6.2f , %6.2f , %6.2f ', s2, S2, T2, c2, tElapsed);
                fprintf(fid,'\n');
            end
        end
    end
end

fclose(fid);
