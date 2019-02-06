function [y S N]=sSTCdetOpt_FixedT(T,Kr,K0,mi,h,p)
a = p./(p+h);
b0 = sqrt(2*K0./(a*h*mi));
if mod(b0./T,1)==0 % b0/T is integer
    N=1;
    S=a.*mi.*T;
    y = sqrt(2.*K0.*h.*a.*mi) + Kr./T;
else
    Nlo = floor(b0./T);
    Slo = Nlo.*T;
    Nhi = ceil(b0./T);
    Shi = Nhi.*T;
    clo = eoq(Slo,K0,a,h,mi);
    chi = eoq(Shi,K0,a,h,mi);
    if clo < chi
        N=Nlo;
        S=Slo;
        y=clo + Kr./T;
    else
        N=Nhi;
        S=Shi;
        y=chi + Kr./T;
    end
end
end

function y=eoq(S,K0,a,h,m)
y = K0./S + a.*h.*m.*S./2;
end
    