function [xb yb xs k] = nonlinoptEA(f, x0, sigma, kmax, quantvec)
if nargin < 5
    quantvec = ones(numel(x0),1);
end
if nargin < 4
    kmax = 1000;
end
n = numel(x0);
xs = ones(n+1,kmax);
x = x0;
xb = x0;
yb = f(x0);
fx = yb;
zerosn = zeros(n,1);
for i=1:kmax
    xtry = x + random('norm',zerosn, sigma);
    for j=1:n
        xtry(j) = round(xtry(j)/quantvec(j,1))*quantvec(j,1);
    end
    ftry = f(xtry);
    if ftry < fx
        x = xtry;
        fx = ftry;
        if ftry < yb
            yb = ftry;
            xb = x;
        end
    end
    xs(1:n,i)=x;
    xs(n+1,i)=fx;
end
k=i;
end

