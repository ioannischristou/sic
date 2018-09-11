% comparison between exact (s,S,T) and approximation algorithms

fid = fopen('sST_minM.txt', 'wt');

L=1;
lamda = 50;
%tnot = 10^-6;
tnot = 5 / lamda;  % require lamda*T >= 5
epst=0.01;
eps=10^-4;
p2 = 0;
estop = 1;
Tmax = 5.0;

% sigma M Kr K0 h p s* S* T* c* time 

fprintf(fid, '     M  ,    Kr  ,     Ko ,      h ,      p ,   s*   ,   S*   ,  T*    ,   cost , time\n');

% Krs must be in descending order!
Krs=[5 2 1 0.1];
Kr_len = numel(Krs);
K0s=[1 5 25 100];
K0_len = numel(K0s);
hs=[10 15 20];
h_len = numel(hs);
ps = [20 25 100];
p_len = numel(ps);
Ms=[10 20 50 100 200];
m_len = numel(Ms);

mi=lamda;
sigma=sqrt(lamda);


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
                tic;
                [s S T c] = sSTCpoissonOpt(3,Kr,K0,L,lamda,h,p,0,0.01,tnot);
                tElapsed=toc;
                % heuristically find a feasible (s,S,T) policy.
                for n=1:m_len
                    Mi=Ms(n);
                    if S-s<Mi
                        S2=s+Mi;
                        c2=sSTCpoisson(s,S2,T,Kr,K0,L,lamda,h,p);
                        s3=S-Mi;
                        c3=sSTCpoisson(s3,S,T,Kr,K0,L,lamda,h,p);
                        if c2<c3
                            S=S2;
                            c=c2;
                        else
                            s=s3;
                            c=c3;
                        end
                    end
                    fprintf(fid,'%6.2f , %6.2f , %6.2f , %6.2f , %6.2f ,', Mi,Kr,K0,h,p);
                    fprintf(fid,'%6.2f , %6.2f , %6.2f , %6.2f , %6.2f \n', s, S, T, c, tElapsed);
                end
            end
        end
    end
end
fclose(fid);
