function [s S T c] = sSTCpoissonOptFastApproxM(Kr,K0,L,lamda,h,p, p2, epst, tdelta, M)
% utilizes the fast approximate search heuristic of LCS for 
% determining the optimal parameter policy for (r,nQ,T) with
% poisson demand (normal demand with mi=lamda, sigma=sqrt(lamda)) and then 
% utilizes the Ti locally optimal times it finds from the
% snQTCnormalOptFastApproxM() method to search for the optimal 
% cost of the (s,S,Ti) policies -using the 
% Zheng & Federgruen algorithm- and picks the best one.
%tic;
if nargin < 7
    p2=0;
end
if nargin < 8
    epst = 0.01;
end
if nargin < 9
    tdelta = 0.1;
end
if nargin < 10
    M=1;
end

mi = lamda;
sigma = sqrt(lamda);
[ri Qi Ti ci] = snQTCnormOptFastApproxM(Kr,K0,L,mi,sigma,h,p,p2,M);
Ti_len = numel(Ti);
for i=1:Ti_len
    disp(['Ti=' num2str(Ti(i))]);
end

copt = 10^30;
for i=1:Ti_len
    T1=Ti(i);
    for T=T1-tdelta:epst:T1+tdelta
        if T<0
            T=epst;
        end
        [s S c] = sSTCpoissonOpt_FixedT(T, Kr, K0, L, lamda, h, p, p2);
        if c < copt
            sopt = s;
            Sopt = S;
            Topt = T;
            copt = c;
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