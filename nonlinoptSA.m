function [xbest ybest xs k] = nonlinoptSA(f, x0, kmax, h, deltagenerator, schedule)
if nargin < 3
    kmax = 1000;
end
if nargin < 6
    schedule = @(x) invt(x, kmax);
end
if nargin < 4
    h = 1;
end
if nargin < 5
    deltagenerator = @(t,h) randdelta(t,h);
end

tnot = 10^-6;

n = numel(x0);
xs = zeros(n+1, kmax);
fx = f(x0);
disp(['starting at f(x0)=' num2str(fx)]);
x = x0;
xbest = x0;
ybest = fx;
acc = 0;
rej = 0;
for i=1:kmax
    T = schedule(i);
    if T < tnot
        break;
    end
    xtry = x + deltagenerator(x,h);
    ftry = f(xtry);
    % keep incumbent
    if ftry < ybest
        xbest = xtry;
        ybest = ftry;
    end
    Df = ftry - fx;
    if Df < 0
        x = xtry;
        fx = ftry;
        acc = acc+1;
    else
        P = random('unif',0,1);
        thres = exp(-Df/T);
        if P <= thres
            x = xtry;
            fx = ftry;
            acc = acc+1;
        else 
            rej = rej+1;
        end
    end
    % keep trajectory
    xs(1:n,i) = x(1:n,1);
    xs(n+1,i) = fx;
end
k = i;
disp(['acc=' num2str(acc) ' rej=' num2str(rej)]);  % itc: HERE rm asap
end

function d = randdelta(x, h)
n = numel(x);
d = zeros(n,1);
for i=1:n
    d(i,1) = random('unif',-h,h);
end
end

function y = invt(k, kmax)
T0 = kmax;  % initial temperature
% a mod d = a-d.*floor(a./d);
lvl = 1+floor(20*k/T0);
y = T0/lvl.^2;
end

