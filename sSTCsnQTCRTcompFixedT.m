% comparison between (R,T) and (r,nQ,T) policies as well as (R,1), (r,nQ,1)
% policies

fid = fopen('sSTsnQTCRT_Rao_T.txt', 'wt');

L=1;
lamda = 50;
%tnot = 10^-6;
tnot = 5 / lamda;  % require lamda*T >= 5
epst=0.01;
p2 = 0;
estop = 1;

% Kr K0 h p s* S* T* c* s* Q* T* c* R* T* c* 

fprintf(fid, 'Kr ,Ko ,h ,p ,s ,S ,c ,r ,Q ,c ,R ,c ,s ,S, c ,r ,Q ,c, R ,c ,s ,S ,c ,r ,Q ,c ,R ,c ,s ,S ,c ,r ,Q ,c ,R ,c \n');

% Krs must be in descending order!
Krs=[0];
Kr_len = numel(Krs);
K0s=[1000];
K0_len = numel(K0s);
hs=[10 15 20];
h_len = numel(hs);
ps = [20 25 100];
p_len = numel(ps);
Ts = [2.0];
t_len = numel(Ts);

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
            for i=1:Kr_len
                Kr=Krs(i);
                disp(['computing (s,nQ,T) and (R,T) policies for Kr=' num2str(Kr)]);
                fprintf(fid,'%6.0f , %6.0f , %6.0f , %6.0f ,', Kr, K0, h, p);
                for m=1:t_len
                    T=Ts(m);
                    [s S c] = sSTCpoissonOpt_FixedT(T,Kr,K0,L,lamda,h,p);
                    fprintf(fid, '%6.2f , %6.2f , %6.2f ,', s, S, c);
                    [sopt Qopt copt] = snQTCpoissonOptFast2_FixedT(T, Kr, K0, L, lamda, h, p);
                    fprintf(fid,'%6.2f , %6.2f , %6.2f , ', sopt, Qopt, copt);
                    [Ropt Q copt] = snQTCpoissonOptFast2_FixedT(T, Kr, K0, L, lamda, h, p, 1, 1);  % run the (R,T) policy
                    fprintf(fid,'%6.2f , %6.2f ,', Ropt, copt);                
                end
                fprintf(fid,'\n');
            end
        end
    end
end

fclose(fid);
