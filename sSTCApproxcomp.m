% comparison between exact (s,S,T) and approximation algorithms

fid = fopen('sSTCApprox.txt', 'wt');

L=1;
lamda = 50;
%tnot = 10^-6;
tnot = 5 / lamda;  % require lamda*T >= 5
epst=0.01;
p2 = 0;
estop = 1;
Tmax = 5.0;

% Kr K0 h p s* S* T* c* time* s1 S1 T1 c1 time1 s2 S2 T2 c2 time2 

fprintf(fid, '    Kr ,     Ko ,      h ,      p ,    s* ,    S* ,    T* ,    c* , time*, s1 , S1 , T1 , c1 , time1, s2 , S2 , T2 , c2 , time2\n');

% Krs must be in descending order!
Krs=[5 2 1 0.1];
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
            [si Si Ti ci timesi] = sSTCpoissonOpt_MultiKr(Tmax, Krs, K0, L, lamda, h, p, p2, epst, tnot, estop);
            % now compute for each Kr the (s,nQ,T) and (R,T) policies
            for i=1:Kr_len
                Kr=Krs(i);
                fprintf(fid,'%6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f ,', Kr, K0, h, p, si(i), Si(i), Ti(i), ci(i), timesi(i));
                disp(['computing Approx policies for M=1 for Kr=' num2str(Kr)]);
                tic;
                [sopt Sopt Topt copt] = sSTCpoissonOptFastApproxM(Kr, K0, L, lamda, h, p, p2, epst, epst, 1);
                tElapsed=toc;
                fprintf(fid,'%6.2f , %6.2f , %6.2f , %6.2f , %6.2f ,', sopt, Sopt, Topt, copt, tElapsed);
                disp(['computing Approx policies for M=2 for Kr=' num2str(Kr)]);
                tic;
                [sopt Sopt Topt copt] = sSTCpoissonOptFastApproxM(Kr, K0, L, lamda, h, p, p2, epst, epst, 2);
                tElapsed=toc;
                fprintf(fid,'%6.2f , %6.2f , %6.2f , %6.2f , %6.2f', sopt, Sopt, Topt, copt, tElapsed);
                fprintf(fid,'\n');
            end
        end
    end
end

fclose(fid);
