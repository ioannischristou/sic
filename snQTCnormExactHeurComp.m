% comparison between exact and heuristic determination of (r,nQ,T) policies 

fid = fopen('snQTCnormComp.new.txt', 'wt');

L=1;
lamda = 50;
mi=lamda;
sigma=sqrt(lamda);
%tnot = 10^-6;
tnot = (3*sigma/mi)^2;
epst=0.01;
epsq=0.1;
qnot=0.001;
p2 = 0;
estop = 1;
Tmax = 5.0;

% Kr K0 h p s* Q* T* c* time s1 Q1 T1 c1 time 

fprintf(fid, '    Kr ,     Ko ,      h ,      p ,    s* ,    Q* ,    T* ,    c* , time  ,     s1 ,     Q1 ,    T1 ,    c1 ,  time \n');

% Krs must be in descending order!
Krs=[5 2 1];
Kr_len = numel(Krs);
K0s=[1 5 25 100 1000];
K0_len = numel(K0s);
hs=[10 15 20];
%hs=[20];
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
            % now compute for each Kr the (s,nQ,T) and (R,T) policies
            for i=1:Kr_len
                Kr=Krs(i);
                disp(['computing exact (s,nQ,T) policies for Kr=' num2str(Kr)]);
                tic;
                [sopt Qopt Topt copt] = snQTCnormOptFast2(Kr, K0, L, mi, sigma, h, p, epsq, epst, qnot, tnot);
                tEl=toc;
                disp(['computing heuristic (s,nQ,T) policies for Kr=' num2str(Kr)]);
                tic;
                [s Q T c1 s2 q2 t2 c2 c] = snQTCnormOptFastApprox2(Kr,K0,L,mi,sigma,h,p);
                tEl2=toc;
                if c < copt  % quantization error
                    sopt=s;
                    Qopt=Q;
                    Topt=T;
                    copt=c;
                end
                fprintf(fid,'%6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f ,', Kr, K0, h, p, sopt, Qopt, Topt, copt, tEl);
                fprintf(fid,'%6.2f , %6.2f , %6.2f , %6.2f , %6.2f', s, Q, T, c, tEl2);
                fprintf(fid,'\n');
            end
        end
    end
end

fclose(fid);
