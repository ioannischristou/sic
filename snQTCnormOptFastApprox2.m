function [sopt qopt topt copt sopt2 qopt2 topt2 copt2 creal] = snQTCnormOptFastApprox2(Kr, K0, L, mi, sigma, h, p, p2, epsq, epst)
% use the approximation for the Porder so as to obtain two convex
% optimization problems and report the min of the two
% the values epsq and epst, if present indicate a quantization of the 
% Q and/or T parameters and are treated in the end of the computation.
if nargin < 9
    epsq = 0;
end
if nargin <10
    epst = 0;
end
if nargin < 8
    p2 = 0;
end
% mi*T-3sigma*sqrt(T)>=0 <=> sqrt(T)>=3sigma/mi
Tmin = (3.5*sigma/mi)^2;

T0 = 1;
if T0<Tmin
    T0=Tmin;
end
s0 = mi*(L+T0);
Q0 = 1;
smin = -10^10;  % used to be -10^30
v0 = [s0, T0];  % initial estimate

[v1, fval1] = fmincon(@(v) f1(v, Kr, K0, L, mi, sigma, h, p, p2), v0, [], [], [], [], [smin, Tmin], []);  % used to be 1.e-8 instead of Tmin
if v1(2)>=100  % hack to get around instability issue with optimization of (R,T) in a system with b(t)=p + 0t
    [v1, fval1] = fmincon(@(v) f11(v, Kr, K0, L, mi, sigma, h, p, p2), v0, [], [], [], [], [smin, Tmin], [-smin, 100]);
end

v0 = [s0, Q0, T0];
[v2, fval2] = fmincon(@(v) f2(v, Kr, K0, L, mi, sigma, h, p, p2), v0, [], [], [], [], [smin, 1.e-9, Tmin], [-smin, 1.e9, 100]);
if fval1 < fval2
    sopt = v1(1);
    qopt = 0;
    topt = v1(2);
    copt = fval1;
    sopt2 = v2(1);
    qopt2 = v2(2);
    topt2 = v2(3);
    copt2 = fval2;
    creal = fval1;
else
    sopt = v2(1);
    qopt = v2(2);
    topt = v2(3);
    copt = fval2;
    sopt2 = v1(1);
    qopt2 = 0;
    topt2 = v1(2);
    copt2 = fval1;
    creal = snQTCnorm2(sopt,qopt,topt,Kr,K0,L,mi,sigma,h,p);
    if p2>0
        creal = creal + p2*((h1(sopt,topt,L,mi,sigma)-h1(sopt+qopt,topt,L,mi,sigma))/qopt);
    end
end
% in case of quantized Q or T, figure out best "nearest" soln
if epsq > 0
    qnew1 = nearest(qopt,epsq);
    qnew2 = nearest(qopt2,epsq);
    if epst > 0
        % both Q&T are quantized
        tnew1 = nearest(topt,epst);
        tnew2 = nearest(topt2, epst);
    else % epst == 0
        % only quantize Q
        tnew1 = topt;
        tnew2 = topt2;
    end
    % let s' be the argmin. of c(s,Q,T)
    s0=mi*(L+tnew1);
    smax = 10^30;
    smin = -10^30;
    [snew1 val]=lsqnonlin(@(x) aTLCeq(x,qnew1,tnew1,L,mi,sigma,h,p), s0, smin, smax);
    if abs(val) > 10^-3
        error(['h=' num2str(h) ' p=' num2str(p) ' Q=' num2str(qnew1) ' T=' num2str(tnew1) ' s=' num2str(snew1) ' not solved val=' num2str(val)]);
    end
    s0 = mi*(L+tnew2);
    [snew2 val]=lsqnonlin(@(x) aTLCeq(x,qnew2,tnew2,L,mi,sigma,h,p), s0, smin, smax);
    if abs(val) > 10^-3
        error(['h=' num2str(h) ' p=' num2str(p) ' Q=' num2str(qnew2) ' T=' num2str(tnew2) ' s=' num2str(snew2) ' not solved val=' num2str(val)]);
    end    
    cnew1 = snQTCnorm2(snew1, qnew1, tnew1, Kr, K0, L, mi, sigma, h, p);
    cnew2 = snQTCnorm2(snew2, qnew2, tnew2, Kr, K0, L, mi, sigma, h, p);
    if cnew1 < cnew2
        sopt = snew1;
        sopt2 = snew2;
        qopt = qnew1;
        qopt2 = qnew2;
        topt = tnew1;
        topt2 = tnew2;
        copt = cnew1;
        copt2 = cnew2;
        creal = copt;
        if p2>0 && qopt>0
           creal = creal + p2*((h1(sopt,topt,L,mi,sigma)-h1(sopt+qopt,topt,L,mi,sigma))/qopt); 
        end
    else
        sopt = snew2;
        sopt2 = snew1;
        qopt = qnew2;
        qopt2 = qnew1;
        topt = tnew2;
        topt2 = tnew1;
        copt = cnew2;
        copt2 = cnew1;
        creal = copt;
        if p2>0 && qopt>0
           creal = creal + p2*((h1(sopt,topt,L,mi,sigma)-h1(sopt+qopt,topt,L,mi,sigma))/qopt); 
        end        
    end
else if epst > 0
        % only quantize T
        qnew1 = qopt;
        qnew2 = qopt2;
        tnew1 = nearest(topt,epst);
        tnew2 = nearest(topt2, epst);
        % let s' be the argmin. of c(s,Q,T)
        s0=mi*(L+T);
        smax = 10^30;
        smin = -10^30;
        [snew1 val]=lsqnonlin(@(x) aTLCeq(x,qnew1,tnew1,L,mi,sigma,h,p), s0, smin, smax);
        if abs(val) > 10^-3
            error(['h=' num2str(h) ' p=' num2str(p) ' Q=' num2str(qnew1) ' T=' num2str(tnew1) ' s=' num2str(snew1) ' not solved val=' num2str(val)]);
        end
        [snew2 val]=lsqnonlin(@(x) aTLCeq(x,qnew2,tnew2,L,mi,sigma,h,p), s0, smin, smax);
        if abs(val) > 10^-3
            error(['h=' num2str(h) ' p=' num2str(p) ' Q=' num2str(qnew2) ' T=' num2str(tnew2) ' s=' num2str(snew2) ' not solved val=' num2str(val)]);
        end    
        cnew1 = snQTCnorm2(snew1, qnew1, tnew1, Kr, K0, L, mi, sigma, h, p);
        cnew2 = snQTCnorm2(snew2, qnew2, tnew2, Kr, K0, L, mi, sigma, h, p);
        if cnew1 < cnew2
            sopt = snew1;
            sopt2 = snew2;
            qopt = qnew1;
            qopt2 = qnew2;
            topt = tnew1;
            topt2 = tnew2;
            copt = cnew1;
            copt2 = cnew2;
            creal = copt;
            if p2>0 && qopt>0
                creal = creal + p2*((h1(sopt,topt,L,mi,sigma)-h1(sopt+qopt,topt,L,mi,sigma))/qopt); 
            end
        else
            sopt = snew2;
            sopt2 = snew1;
            qopt = qnew2;
            qopt2 = qnew1;
            topt = tnew2;
            topt2 = tnew1;
            copt = cnew2;
            copt2 = cnew1;
            creal = copt;
            if p2>0 && qopt>0
                creal = creal + p2*((h1(sopt,topt,L,mi,sigma)-h1(sopt+qopt,topt,L,mi,sigma))/qopt); 
            end        
        end    
    end
end

end  % function snQTCnormOptFastApprox2


function y = nearest(x,eps)
if eps==0
    y=x;
    return;
end
d = floor(x/eps);
y = eps*d;
y2 = eps*(d+1);
if abs(y2-x)<abs(y-x)
    y = y2;
end
if y==0
    y = eps;
end
end


function y = f1(v, Kr, K0, L, mi, sigma, h, p, p2)
y = STCnormVec(v, Kr, K0, L, mi, sigma, h, p, p2);
end

function z = f2(v, Kr, K0, L, mi, sigma, h, p, p2)

s = v(1);
Q = v(2);
T = v(3);

P0 = mi*T./Q;

Pi1 = h*(s+Q/2-mi*(L+T/2));

Pi2 = quadgk(@(t) penint(t,s,Q,L,mi,sigma), 0, T);
Pi2 = ((h+p)/(2*T))*Pi2; 
Pi3 = p2*((h1(s,T,L,mi,sigma)-h1(s+Q,T,L,mi,sigma))/Q);
z = (Kr+K0*P0)/T + Pi1 + Pi2 + Pi3;
end

function z = f11(v, Kr, K0, L, mi, sigma, h, p, p2)
s = v(1);
Q = 1.e-6;
T = v(2);

P0 = 1;

Pi1 = h*(s+Q/2-mi*(L+T/2));

Pi2 = quadgk(@(t) penint(t,s,Q,L,mi,sigma), 0, T);
Pi2 = ((h+p)/(2*T))*Pi2; 
Pi3 = p2*((h1(s,T,L,mi,sigma)-h1(s+Q,T,L,mi,sigma))/Q);
z = (Kr+K0*P0)/T + Pi1 + Pi2 + Pi3;
end


function z = penint(T,s,Q,L,mi,sigma)
ZLpt = (s-mi*(L+T))./(sigma*sqrt(L+T));
RLpt = Q./(sigma*sqrt(L+T));
t1=((ZLpt+RLpt).^2+1).*normcdf(ZLpt+RLpt)+(ZLpt+RLpt).*normpdf(ZLpt+RLpt);
t2=(ZLpt.^2+1).*normcdf(ZLpt)+ZLpt.*normpdf(ZLpt)+(2*ZLpt+RLpt).*RLpt;
z = (t1 - t2).*sigma.*sqrt(L+T)./RLpt; 
end


function z=h1(r,T,L,mi,sigma)
z = (l1(r,L+T,mi,sigma)-l1(r,L,mi,sigma))/T;
end

function z = l1(v,T,lamda,D)
y=(v-lamda*T)/sqrt(D*T);
z = ((1+(v-lamda*T)^2/(D*T))*(1-normcdf(y))-y*normpdf(y))*D*T/2;
end