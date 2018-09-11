function [s S T c] = sSTCpoissonOptFastApprox3(Kr,K0,L,lamda,h,p, p2, ddelta, epst, tdelta)
% utilizes the fast approximate search heuristic of LCS for 
% determining the optimal parameter policy for (r,nQ,T) with
% poisson demand and then 
% utilizes the T1 and T2 locally optimal times it finds from the
% snQTCpoissonOptFastApprox() method to search for the optimal 
% cost of the (s,S,T1) and (s,S,T2) policies -using the 
% Zheng & Federgruen algorithm- and picks the best one.
tic;
if nargin < 7
    p2=0;
end
if nargin < 8
    ddelta = 2;
end

if nargin < 9
    epst = 0.01;
end
if nargin < 10
    tdelta = 0.02;
end

[r1 Q1 T1 c_1 r2 Q2 T2 c_2] = snQTCpoissonOptFastApprox3(Kr,K0,L,lamda,h,p,p2);
disp(['T1=' num2str(T1) ' T2=' num2str(T2)]);

c1 = 10^30;
s0 = r1;
S0 = r1+Q1;
T0 = T1;
minT = 5/lamda;
if minT<epst
    minT = epst;
end
% search around s1,S1,T1
for T=T0-tdelta:epst:T0+tdelta
    if T<0
        T=minT;
    end
    %[s S c] = sSTCpoissonOpt_FixedT(T, Kr, K0, L, lamda, h, p, p2);
    for s=s0-ddelta:s0+ddelta
        for S=S0-ddelta:S0+ddelta
            if S<=s
                continue;
            end
            c = sSTCpoisson(s,S,T,Kr,K0,L,lamda,h,p,p2);
            if c < c1
                s1 = s;
                S1 = S;
                T1 = T;
                c1 = c;
            end
        end
    end
end

c2 = 10^30;
% search around T2
s0 = r2;
S0 = r2+Q2;
T0 = T2;
for T=T0-tdelta:epst:T0+tdelta
    if T<0
        T=minT;
    end
    %[s S c] = sSTCpoissonOpt_FixedT(T, Kr, K0, L, lamda, h, p, p2);
    for s=s0-ddelta:s0+ddelta
        for S=S0-ddelta:S0+ddelta
            if S<=s
                continue;
            end
            c = sSTCpoisson(s,S,T,Kr,K0,L,lamda,h,p,p2);
            if c < c2
                s2 = s;
                S2 = S;
                T2 = T;
                c2 = c;
            end
        end
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