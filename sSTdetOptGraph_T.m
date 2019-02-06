function yi = sSTdetOptGraph_T(tmin,tmax,Kr,K0,mi,h,p,tstep)
Ti = tmin:tstep:tmax;
len = numel(Ti);
yi = 1:len;
for i=1:len
    T=Ti(i);
    yi(i)= sSTCdetOpt_FixedT(T,Kr,K0,mi,h,p);
end

hold on
plot(Ti,yi,'k-');
hold off
