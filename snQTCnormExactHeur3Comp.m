% comparison between exact and heuristic determination of (r,nQ,T) policies 

fid = fopen('snQTCnormComp3.m10.txt', 'wt');

L=1;
lamda = 10;
mi=lamda;
sigma=sqrt(lamda);
%tnot = 10^-6;
tnot = (3.5*sigma/mi)^2;
epst=0.01;
epsq=0.1;
eps_s=0.001;
qnot=10^-8;
p2 = 0;
estop = 1;
Tmax = 100.0;

% Kr K0 h p s* Q* T* c* time s1 Q1 T1 c1 time 

fprintf(fid, '    Kr ,     Ko ,      h ,      p ,    s* ,    Q* ,    T* ,    c* , time  ,     s1 ,     Q1 ,    T1 ,    c1 ,  time \n');

% Krs must be in descending order!
Krs=[3 2 0];
Kr_len = numel(Krs);
K0s = [5 25 100 1000];
K0_len = numel(K0s);
hs = [10 15 20];
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
                disp(['computing exact (s,nQ,T) policies for Kr=' num2str(Kr) ' K0=' num2str(K0) ' h=' num2str(h) ' p=' num2str(p)]);
                tic;
                [sopt Qopt Topt copt Qmax Tmax copt2] = snQTCnormOptFast2(Kr, K0, L, mi, sigma, h, p, epsq, epst, qnot, tnot, eps_s);
                tEl=toc;
                disp(['computing heuristic (s,nQ,T) policies for Kr=' num2str(Kr) ' K0=' num2str(K0) ' h=' num2str(h) ' p=' num2str(p)]);
                tic;
                [s Q T c] = snQTCnormOptFastApprox3(Kr,K0,L,mi,sigma,h,p, epsq, epst, qnot, tnot, eps_s);
                tEl2=toc;
                if c < copt  % quantization error
                    disp(['copt=' num2str(copt) ' > capprox=' num2str(c) ' Qapprox=' num2str(Q) ' sapprox=' num2str(s) ' Tapprox=' num2str(T) '|| qopt=' num2str(Qopt) ' sopt=' num2str(sopt) ' Topt=' num2str(Topt)]);
                    disp(['exact method: Qmax=' num2str(Qmax) ' Tmax=' num2str(Tmax)]);
                    if 100*(copt-c)/c > 5.e-2   %c + 1.e-2 < copt  % error is too big!
                        error('must correct this');
                    end
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
