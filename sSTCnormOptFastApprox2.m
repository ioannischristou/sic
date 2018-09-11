function [s S T c] = sSTCnormOptFastApprox2(Kr,K0,L,mi,sigma,h,p, p2, epsd, tdelta)
% utilizes the fast approximate search heuristic of LCS for 
% determining the optimal parameter policy for (r,nQ,T) with
% normal demand (normal demand N(mi*t,t*sigma^2)) and then 
% searches in a neighborhood of the two solutions using standard
% constrained NLP in the form of the fmincon() function.
tic;
if nargin < 8
    p2=0;
end
if nargin < 9
    epsd = 5;
end
if nargin < 10
    tdelta = 0.05;
end

Tmin = (3*sigma/mi)^2;
Dtol = 10^-3;

[r1 Q1 T1 c1 r2 Q2 T2 c2] = snQTCnormOptFastApprox2(Kr,K0,L,mi,sigma,h,p,p2);
if Q1<=0
    Q1 = Dtol;
end
if Q2<=0
    Q2 = Dtol;
end
disp(['T1=' num2str(T1) ' T2=' num2str(T2)]);

A = [1, -1, 0];
b = -Dtol;

s1 = r1;
S1 = r1+Q1;
% search around (s1, S1, T1)
v0 = [s1; S1; T1];
Tlow = T1-tdelta;
Tupp = T1+tdelta;
if Tlow<Tmin
    Tlow = Tmin;
end
if Tupp<Tlow
    Tupp = Tlow+tdelta;
end
lb = [s1-epsd; S1-epsd; Tlow];
ub = [s1+epsd; S1+epsd; Tupp];
[v1 c1] = fmincon(@(v) fP(v,Kr,K0,L,mi,sigma,h,p,p2), v0, A, b, [], [], lb, ub);

s2 = r2;
S2 = r2+Q2;
% search around (s2, S2, T2)
v0 = [s2; S2; T2];
Tlow = T2-tdelta;
Tupp = T2+tdelta;
if Tlow<Tmin
    Tlow = Tmin;
end
if Tupp<Tlow
    Tupp = Tlow+tdelta;
end
lb = [s2-epsd; S2-epsd; Tlow];
ub = [s2+epsd; S2+epsd; Tupp];
[v2 c2] = fmincon(@(v) fP(v,Kr,K0,L,mi,sigma,h,p,p2), v0, A, b, [], [], lb, ub);

if c1 < c2
    s = v1(1);
    S = v1(2);
    T = v1(3);
    c = c1;
else
    s = v2(1);
    S = v2(2);
    T = v2(3);
    c = c2;
end

tElapsed=toc;
disp(['run-time=' num2str(tElapsed) ' (secs)']);
end

function y=fP(v, Kr, K0, L, mi, sigma, h, p, p2)
s = v(1,1);
S = v(2,1);
T = v(3,1);
y=sSTCnorm(s,S,T,Kr,K0,L,mi,sigma,h,p,p2);
%disp(['s=' num2str(s) ' S=' num2str(S) ' T=' num2str(T) ' c=' num2str(y)]);
end

