function y = ordercost(Qmax,K0,T,mi,sigma,epsq,qnot)
% plot order costs for the normal demand case of (r,nQ,T)
if nargin < 7
    qnot = 1.e-8;
end
if nargin < 6
    epsq = 0.1;
end
Qlen = Qmax/epsq+1;

Qi = 1:Qlen;
Ci = 1:Qlen;
CLi = 1:Qlen;
Capproxi = 1:Qlen;

diffapproxi = 1:Qlen;

for i=1:Qlen
    if i==1
        Q=qnot;
    else
        Q=(i-1)*epsq;
    end
    
    Qi(i)=Q;
    
    if Q<=1.e-6
        P0=1;
        PL=1;
        Pa=1;
    else
          qltsrt = (Q-mi.*T)./(sigma.*sqrt(T));
          npdf = normpdf(qltsrt);
          ncdf = normcdf(qltsrt);
          P0 = mi.*T.*ncdf./Q + (1-ncdf) - sigma.*sqrt(T).*npdf./Q;
          % Lagodimos formula for P0
          RT=Q/(sigma*sqrt(T));
          MT=T*mi/(sigma*sqrt(T));
          PL=1-(1/RT)*(normpdf(RT-MT)+(RT-MT)*normcdf(RT-MT)-(normpdf(-MT)-MT*normcdf(-MT)));
          Pa=min(mi.*T./Q,1);
    end
    Ci(i) = K0*P0/T;
    CLi(i) = K0*PL/T;
    Capproxi(i) = K0*Pa/T;
    diffapproxi(i) = Capproxi(i)-CLi(i);
    if diffapproxi(i)<0  % numerical inaccuracies may make this value like -1.e-12
        diffapproxi(i)=0;
    end
end

subplot(2,2,1);
hold on
plot(Qi,Capproxi,'k-');
hold off
hold on
plot(Qi,CLi,'b-');
hold off
subplot(2,2,3);
hold off
plot(Qi,diffapproxi,'r-');
hold off

% now try T+1 on the right column of the plot figure

T=T+1;
for i=1:Qlen
    if i==1
        Q=qnot;
    else
        Q=(i-1)*epsq;
    end
    Qi(i)=Q;
    if Q<=1.e-6
        P0=1;
        PL=1;
        Pa=1;
    else
          qltsrt = (Q-mi.*T)./(sigma.*sqrt(T));
          npdf = normpdf(qltsrt);
          ncdf = normcdf(qltsrt);
          P0 = mi.*T.*ncdf./Q + (1-ncdf) - sigma.*sqrt(T).*npdf./Q;
          % Lagodimos formula for P0
          RT=Q/(sigma*sqrt(T));
          MT=T*mi/(sigma*sqrt(T));
          PL=1-(1/RT)*(normpdf(RT-MT)+(RT-MT)*normcdf(RT-MT)-(normpdf(-MT)-MT*normcdf(-MT)));
          Pa=min(mi.*T./Q,1);
    end
    Ci(i) = P0;
    CLi(i) = PL;
    Capproxi(i) = Pa;
    diffapproxi(i) = Capproxi(i)-CLi(i);
    
end

subplot(2,2,2);
hold on
plot(Qi,Capproxi,'k-');
hold off
hold on
plot(Qi,CLi,'b-');
hold off
subplot(2,2,4);
hold off
plot(Qi,diffapproxi,'r-');
hold off


end
