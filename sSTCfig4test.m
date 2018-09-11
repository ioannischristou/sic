function y=sSTCfig4test()
rp=3;
%figure('Fig4');
for j=1:2*rp, subplot(rp,2,j)
    % figure out the three graphs
    [x y] = STgraph(j);
    plot(x,y,'b');
    [x y] = snQTgraph(j);
    hold on;
    plot(x,y,'r');
    [x y] = sSTgraph(j);
    hold on;
    plot(x,y,'k');
end
end

function [x y]=STgraph(j)
L=0;
h=1;
p=9;
% figure out K0 and lamda
if j<=2
    K0=64;
    l=10;
else if j<=4
        K0=32;
        l=5;    
    else
        K0=64;
        l=46;
    end
end
% figure out Kr
if j==1 || j==3 || j==5
    Kr=0.1;
    epst=0.05;
    Tmax=2.0;
else
    Kr=K0/2;
    epst=0.1;
    Tmax=7;
end
% now get the cost for T=0.1:0.05:2
Ti=0.1:epst:Tmax;
Tlen = numel(Ti);
[ds dq dt costs] = snQTCpoissonOptGraph_T(1,Tmax,Kr,K0,L,l,h,p,1,epst);
y = 1:Tlen;
for i=1:Tlen
    y(i)=costs(1,i);
end
x=Ti;
end

function [x y] = snQTgraph(j)
L=0;
h=1;
p=9;
% figure out K0 and lamda
if j<=2
    K0=64;
    l=10;
else if j<=4
        K0=32;
        l=5;    
    else
        K0=64;
        l=46;
    end
end
% figure out Kr
if j==1 || j==3 || j==5
    Kr=0;
    epst=0.05;
    Tmax=2.0;
else
    Kr=K0/2;
    epst=0.1;
    Tmax=7;
end
% now get the cost for T=0.1:0.05:2
%Ti=0.1:epst:Tmax;
[ds dq dt x y] = snQTCpoissonOptGraph_T2(Tmax,Kr,K0,L,l,h,p,1,epst,0.1);
end

function [x y]=sSTgraph(j)
L=0;
h=1;
p=9;
% figure out K0 and lamda
if j<=2
    K0=64;
    l=10;
else if j<=4
        K0=32;
        l=5;    
    else
        K0=64;
        l=46;
    end
end
% figure out Kr
if j==1 || j==3 || j==5
    Kr=0;
    epst=0.05;
    Tmax=2;
else
    Kr=K0/2;
    epst=0.1;
    Tmax=7;
end
% now get the cost for T=0.1:0.05:2
[ds dS dT dc dsi dSi dfl x y] = sSTCpoissonOpt(Tmax,Kr,K0,L,l,h,p,0,epst,0.1,0);

end

