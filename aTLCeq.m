function y=aTLCeq(s,Q,T,L,mi,sigma,h,p)

y = quad(@(t) aTLC(t,s,Q,L,mi,sigma),0,T);
y = y/T;

y = y - p/(p+h);
end

function z = aTLC(t,s,Q,L,mi,sigma)
RLpt = Q./(sigma.*sqrt(L+t));
ZLpt = (s-mi*(L+t))./(sigma.*sqrt(L+t));
z = (normpdf(ZLpt+RLpt) + (ZLpt+RLpt).*normcdf(ZLpt+RLpt) - normpdf(ZLpt) - ZLpt.*normcdf(ZLpt))./RLpt;
end

