function y=snQTCnorm2(s,Q,T,Kr,K0,L,mi,sigma,h,p)

RT=Q/(sigma*sqrt(T));
MT=T*mi/(sigma*sqrt(T));
P0=1-(1/RT)*(normpdf(RT-MT)+(RT-MT)*normcdf(RT-MT)-(normpdf(-MT)-MT*normcdf(-MT)));
% Pi1=h*(s+Q/2-mi*(L+(T+1)/2));
Pi1 = h*(s+Q/2-mi*(L+T/2));

Pi2 = quadgk(@(t) penint(t,s,Q,L,mi,sigma), 0, T);  % was quad
Pi2 = ((h+p)/(2*T))*Pi2;  
y = (Kr+K0*P0)/T + Pi1 + Pi2;
end


function z = penint(T,s,Q,L,mi,sigma)
ZLpt = (s-mi*(L+T))./(sigma*sqrt(L+T));
RLpt = Q./(sigma*sqrt(L+T));
t1=((ZLpt+RLpt).^2+1).*normcdf(ZLpt+RLpt)+(ZLpt+RLpt).*normpdf(ZLpt+RLpt);
t2=(ZLpt.^2+1).*normcdf(ZLpt)+ZLpt.*normpdf(ZLpt)+(2*ZLpt+RLpt).*RLpt;
z = (t1 - t2).*sigma.*sqrt(L+T)./RLpt; 
end

