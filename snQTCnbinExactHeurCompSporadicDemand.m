% comparison between exact and heuristic determination of (r,nQ,T) policies 
% itc20180910: file name used to be snQTCnbinCompSporadicDemand.new.txt
fid = fopen('snQTCnbinCompSporadicDemand.new2.txt', 'wt');

L=1;
%lamdas=[1 5 10 15 20 25 26 27 28 29 30 35 40];
lamdas=[0.5];
pl=0.8;  % large variance
epst=0.01;
epsq=1;
estop = 0;
tmin = 0.0001;
mint = 0;
minq = 0;

% Kr K0 h p s* Q* T* c* time s1 Q1 T1 c1 time 

fprintf(fid, '    Kr ,     Ko ,  lamda ,    h ,      p ,    s* ,    Q* ,    T* ,    c* , time  ,     s1 ,     Q1 ,    T1 ,    c1 ,  time \n');

% Krs must be in descending order!
Krs=[1 0.1 0.05];
Kr_len = numel(Krs);
K0s=[1 5 25 100];
K0_len = numel(K0s);
hs=[10 15 20];
%hs=[1];
h_len = numel(hs);
%ps = [1 9];
ps = [20 25 100];
p_len = numel(ps);
lamda_len = numel(lamdas);
p2 = 0;

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
                for m=1:lamda_len
                    lamda=lamdas(m);
                    disp(['computing exact (s,nQ,T) policies for m=' num2str(m)]);
                    tic;
                    [sopt Qopt Topt copt] = snQTCnbinOptFast3SD(minq, mint, Kr, K0, L, lamda, pl, h, p, epsq, epst, estop, tmin);
                    tEl=toc;
                    disp(['computing heuristic (s,nQ,T) policies for m=' num2str(m)]);
                    tic;
                    [s Q T c] = snQTCnbinOptFastApprox5(Kr,K0,L,lamda,pl,h,p,p2,tmin);
                    tEl2=toc;
                    if c < copt  % quantization error
                        sopt=s;
                        Qopt=Q;
                        Topt=T;
                        copt=c;
                    end
                    fprintf(fid,'%6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f , %6.2f ,', Kr, K0, lamda, h, p, sopt, Qopt, Topt, copt, tEl);
                    fprintf(fid,'%6.2f , %6.2f , %6.2f , %6.2f , %6.2f', s, Q, T, c, tEl2);
                    fprintf(fid,'\n');
                end
            end
        end
    end
end

fclose(fid);
