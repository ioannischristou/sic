function [s S T c] = sSTCnormOptFastApprox(Kr,K0,L,mi,sigma,h,p, p2, epsd, epst, tdelta)
% utilizes the fast approximate search heuristic of LCS for 
% determining the optimal parameter policy for (r,nQ,T) with
% normal demand (normal demand N(mi*t,t*sigma^2)) and then 
% utilizes the T1 and T2 locally optimal times it finds from the
% snQTCnormalOptFastApprox() method to search for the optimal 
% cost of the (s,S,T1) and (s,S,T2) policies -using the 
% Zheng & Federgruen algorithm- and picks the best one.
tic;
if nargin < 8
    p2=0;
end
if nargin < 9
    epsd = 1;
end
if nargin < 10
    epst = 0.01;
end
if nargin < 11
    tdelta = 0.1;
end

[r1 Q1 T1 c1 r2 Q2 T2 c2] = snQTCnormOptFastApprox2(Kr,K0,L,mi,sigma,h,p,p2);
disp(['T1=' num2str(T1) ' T2=' num2str(T2)]);

c1 = 10^30;
% search around T1
for T=T1-tdelta:epst:T1+tdelta
    if T<0
        T=epst;
    end
    [s S c] = sSTCnormOpt_FixedT(T, Kr, K0, L, mi, sigma, h, p, p2, epsd);
    if c < c1
        s1 = s;
        S1 = S;
        T1 = T;
        c1 = c;
    end
end

c2 = 10^30;
% search around T2
for T=T2-tdelta:epst:T2+tdelta
    if T<0
        T=epst;
    end
    [s S c] = sSTCnormOpt_FixedT(T, Kr, K0, L, mi, sigma, h, p, p2, epsd);
    if c < c2
        s2 = s;
        S2 = S;
        T2 = T;
        c2 = c;
    end
end

if c1 < c2
    s = s1;
    S = S1;
    T = T1;
    c = c1;
else
    s = s2;
    S = S2;
    T = T2;
    c = c2;
end

tElapsed=toc;
disp(['run-time=' num2str(tElapsed) ' (secs)']);
end