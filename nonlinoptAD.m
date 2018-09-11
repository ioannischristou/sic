function [xopt yopt k] = nonlinoptAD(f, x0, tryorder, epsx, kmax)
% compute a local minimum of the function f starting from x0,
% by minimizing the function f along one of its variables every time
if nargin < 5
    kmax = 10^10;
end
n = numel(x0);
if nargin < 4
    epsx = 1:n;
    for i=1:n
        epsx(i)=0;
    end
end
if nargin < 3
    tryorder = 1:n;
end

x = x0;
k=0;
xopt = x0;
yopt = f(xopt);
yb = yopt;
nto = numel(tryorder);
while k<kmax
    cont = false;
    for i=1:nto
        j=tryorder(i);
        epsxj = epsx(j);
        [xj yj] = optimize1D(f, x, j, epsxj);
        k = k+1;
        if yj<yb
            disp(['nonlinoptAD(): found better solution w/ x(' num2str(j) ')=' num2str(xj) ' -> ' num2str(yj)]);
            cont = true;
            x(j)=xj;
            yb = yj;
        end
    end
    if cont==false
        break;
    end
end
xopt = x;
yopt = yb;
disp(['nonlinoptAD(): done']);
end


function y = f_one(f, x0, j, xj)
x = x0;
x(j) = xj;
y = f(x);
end


function [x y] = optimize1D(f, x, j, epsxj)
myf = @(xj) f_one(f, x, j, xj);
if abs(epsxj)<10^-6
    [x y] = fminunc(myf,x(j));
    return;
end
% find first minimizer of f along x_j direction taking steps of length
% epsxj.
x = findminx(myf, x(j), epsxj);
y = myf(x);
end


function sqt = findminx(f_one, s0, epsx)
step = epsx;  % always start with step-size epsx
niterbnd = 5;
mul = 2;
s = s0;
sqt = s;
cnt = niterbnd;
prevdir=0;
while true
    cnt = cnt-1;
    if cnt==0
        step=step*mul;
        %disp(['step=' num2str(step)]);
        cnt=niterbnd;
    end
    [dir news] = detdir(f_one, s, epsx);
    %if news ~= s
    %    disp(['s=' num2str(s) ' news=' num2str(news)]);  % itc: HERE rm asap
    %end
    if dir==0
        sqt = news;
        break;
    end
    if dir>0
        if prevdir<0 && step>epsx
            step = step/mul;
            cnt=niterbnd;
            %disp(['step=' num2str(step)]);
        end
        s = s+step;  % itc: HERE this should probably go before if stmt
    else % dir<0
        s = s-step;
        if prevdir>0 && step>epsx
            step = step/mul;
            cnt=niterbnd;
            %disp(['step=' num2str(step)]);
        end
    end
    %disp(['s=' num2str(s) ' dir=' num2str(dir) ' prevdir=' num2str(prevdir)]);    
    prevdir=dir;
end
end

function [dir news] = detdir(f_one, s, eps)
c = f_one(s);
cup=f_one(s+eps);
if c > cup
    dir=1;
    news=s;
    return;
else if c==cup
        s2=s;
        while f_one(s2)==c
            %disp(['detdir: up(+)s2=' num2str(s2) ' c=' num2str(c)]);  % itc: HERE rm asap
            s2=s2+eps;
        end
        cnew = f_one(s2);
        if cnew < c
            dir=1;
            news=s2;
            return;
        else % cnew > c
            s2 = s;
            while f_one(s2)==c
                %disp(['detdir: up(-)s2=' num2str(s2) ' c=' num2str(c)]);  % itc: HERE rm asap
                s2 = s2-eps;
            end
            cnew = f_one(s2);
            if cnew < c
                dir = -1;
                news=s2;
                return;
            else % cnew > c
                dir = 0;
                news = s;
                return;
            end
        end
    end
end
cdown = f_one(s-eps);
if c > cdown
    dir=-1;
    news=s;
    return;
else if c==cdown
        s2=s;
        while f_one(s2)==c
            %disp(['detdir: down(-)s2=' num2str(s2) ' c=' num2str(c)]);  % itc: HERE rm asap
            s2=s2-eps;
        end
        cnew = f_one(s2);
        if cnew < c
            dir=-1;
            news=s2;
            return;
        else % cnew > c
            s2 = s;
            while f_one(s2)==c
                %disp(['detdir: down(+)s2=' num2str(s2) ' c=' num2str(c)]);  % itc: HERE rm asap
                s2 = s2+eps;
            end
            cnew = f_one(s2);
            if cnew < c
                dir = 1;
                news=s2;
                return;
            else % cnew > c
                dir = 0;
                news = s;
                return;
            end
        end
    end
end
% if we reach here we're optimal
dir=0;
news=s;
end

