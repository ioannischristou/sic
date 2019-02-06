function [si qi ti] = snQTCnormOptGraph_sigma(sigmamax,Kr,K0,L,mi,h,p, epssigma)
if nargin<8
    epssigma=1;
end
slen=sigmamax/epssigma;
xi = 1:slen;
si = 1:slen;
qi=1:slen;
ti=1:slen;
i=1;
for sigma=1:epssigma:sigmamax
    [sopt qopt topt copt] = snQTCnormOptFastApprox2(Kr,K0,L,mi,sigma,h,p);
    si(i)=sopt;
    qi(i)=qopt;
    ti(i)=topt;
    xi(i)=sigma;
    i = i+1;
end

% plot data
hold on
plot(xi,si,'color',[0.8 0.8 0.8]);
hold on
plot(xi,qi,'b');
plot(xi,ti,'k');

return

