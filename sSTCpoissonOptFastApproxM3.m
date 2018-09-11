function [s S T c] = sSTCpoissonOptFastApproxM3(Kr,K0,L,lamda,h,p, p2, epsd, epst, tdelta, M)
% utilizes the fast approximate search heuristic of LCS for 
% determining the optimal parameter policy for (r,nQ,T) with
% poisson demand and then 
% utilizes the Ti locally optimal times it finds from the
% snQTCpoissonOptFastApproxM3() method to search for the optimal 
% cost of the (s,S,Ti) policies and picks the best one.
%tic;
if nargin < 7
    p2=0;
end
if nargin < 8
    epsd = 5;
end
if nargin < 9
    epst = 0.01;
end
if nargin < 10
    tdelta = 0.1;
end
if nargin < 11
    M=1;
end

[ri Qi Ti ci] = snQTCpoissonOptFastApproxM(Kr,K0,L,lamda,h,p,p2,M);
Ti_len = numel(Ti);
for i=1:Ti_len
    disp(['Ti=' num2str(Ti(i))]);
end

copt = 10^30;
for i=1:Ti_len
    T1=Ti(i);
    s1 = ri(i);
    S1 = ri(i)+Qi(i);
    disp(['start with: s1=' num2str(s1) ', S1=' num2str(S1) ', T1=' num2str(T1)]);
    for T=T1-tdelta:epst:T1+tdelta
        if T<0
            T=epst;
        end
        %[s S c] = sSTCpoissonOpt_FixedT(T, Kr, K0, L, lamda, h, p, p2);
        for s=s1:s1+epsd
            for S=S1:S1+epsd
                if S<=s
                    S=s+1;
                end
                disp(['computing sSTCpoisson( ' num2str(s) ' , ' num2str(S) ' , ' num2str(T) ' )']);
                c = sSTCpoisson(s,S,T,Kr,K0,L,lamda,h,p,p2);
                if c < copt
                    sopt = s;
                    Sopt = S;
                    Topt = T;
                    copt = c;
                end
            end
        end
    end
end
s=sopt;
S=Sopt;
T=Topt;
c=copt;
%tElapsed=toc;
%disp(['run-time=' num2str(tElapsed) ' (secs)']);
end