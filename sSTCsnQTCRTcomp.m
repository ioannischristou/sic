% comparison between (R,T), (r,nQ,T) and (s,S,T)
% policies

fid = fopen('sSTsnQTCRT_Rao.txt', 'wt');

L=1;
lamda = 50;
%tnot = 10^-6;
tnot = 5 / lamda;  % require lamda*T >= 5
epst=0.01;
p2 = 0;
estop = 1;
Tmax = 5.0;

% Kr K0 h p s* S* T* c* s* Q* T* c* R* T* c* 

fprintf(fid, '    Kr ,     Ko ,      h ,      p ,    s* ,    S* ,    T* ,    c* ,     s* ,     Q* ,     T* ,    c* ,    R* ,    T* ,    c*  \n');

% Krs must be in descending order!
Krs=[10 5 3 2 1];
Kr_len = numel(Krs);
K0s=[1000];
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
            [si Si Ti ci] = sSTCpoissonOpt_MultiKr(Tmax, Krs, K0, L, lamda, h, p, p2, epst, tnot, estop);
            % now compute for each Kr the (s,nQ,T) and (R,T) policies
            for i=1:Kr_len
                Kr=Krs(i);
                disp(['computing (s,nQ,T) and (R,T) policies for Kr=' num2str(Kr)]);
                fprintf(fid,'%6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f ,', Kr, K0, h, p, si(i), Si(i), Ti(i), ci(i));
                [sopt Qopt Topt copt] = snQTCpoissonOptFast2(Kr, K0, L, lamda, h, p, 1, epst);
                fprintf(fid,'%6.2f , %6.2f , %6.2f , %6.2f ,', sopt, Qopt, Topt, copt);
                [Ropt Q Topt copt] = snQTCpoissonOptFast2(Kr, K0, L, lamda, h, p, 1, epst, 1);  % run the (R,T) policy
                fprintf(fid,'%6.2f , %6.2f , %6.2f ', Ropt, Topt, copt);                
                fprintf(fid,'\n');
            end
        end
    end
end

fclose(fid);
