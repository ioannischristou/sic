function y=finner3plot(r,R,T,m,s,n,stepsize)
xi=0:stepsize:(R-r);
yi=1:numel(xi);
zi=1:numel(xi);
for i=1:numel(xi)
    x = xi(i);
    yi(i)=finner3(x,r,R,T,m,s,n);
    if x>0
        zi(i)=quadvec(@(u) finner3(u,r,R,T,m,s,n), 0, x);
        if zi(i)<10^-4
            %break up intervals
           num_ints = 50;
           sz = x./num_ints;
           ll=0;
           ul=sz;
           sum=0;
           for j=1:num_ints
               saux = quadvec(@(u) finner3(u,r,R,T,m,s,n), ll, ul);
                if ~isnan(saux)
                    sum = sum + saux;
                end
                ll=ul;
                ul = ul+sz;
           end
           zi(i)=sum;
        end
    else
        zi(i)=0;
    end
end
hold on
plot(xi,yi,'k-');
hold off
hold on
plot(xi,zi,'b-');
end

function y = finner3(x,r,R,T,l,sigma,n)
y1 = n.*complnormcdf(x, l.*T, sigma.*sqrt(T));
y2 = normpdf(R-r-x, (n-1).*l.*T, sigma.*sqrt((n-1).*T));
y = y1*y2;
end

function y = complnormcdf(x,m,s)
y=1.-normcdf(x,m,s);
end
