% for a number of different sigma's, compute the optimal (s,S,T) policy as
% a function of T and plot them all against T.

sigma4 = 4;
sigma3 = 3;
sigma1 = 1; 
sigma01 = 0.1;
Kr=0; % used to be 10
K0=64; % used to be 64
mi=10;
L=0;
h=1; 
p=9; 

%tmin = (sigma3*3/mi)^2;
tmin = 0.01;
tmax = 6.0;
%tstep = 0.025;
tstep = 0.01;

Ti=tmin:tstep:tmax;
Ti_len = numel(Ti);
c3i = 1:Ti_len;
c1i = 1:Ti_len;
c01i = 1:Ti_len;
cpi = 1:Ti_len;
cnbi = 1:Ti_len;

% compute r,p for Nbin given mi=10, s=2
p_nbin = 1-mi./sigma4.^2;
l_nbin = -(1-p_nbin).*log(1-p_nbin).*mi./p_nbin;

% for i=1:Ti_len
%     T = Ti(i);
%     disp(['working on T=' num2str(T)]);
%     [s S c3i(i)] = sSTCnormOpt_FixedT(T,Kr,K0,L,mi,sigma3,h,p,0,0.2);
%     [s S c1i(i)] = sSTCnormOpt_FixedT(T,Kr,K0,L,mi,sigma1,h,p,0,0.2);
%     [s S c01i(i)] = sSTCnormOpt_FixedT(T,Kr,K0,L,mi,sigma01,h,p,0,0.2);
%     [s S cpi(i)] = sSTCpoissonOpt_FixedT(T,Kr,K0,L,mi,h,p);
%     [s Q cnbi(i)] = snQTCnbinOptFast_FixedT(T,Kr,K0,L,l_nbin,p_nbin,h,p);
% end

% hold on
% plot(Ti,c3i,'r-');
% hold off
% hold on
% plot(Ti,c1i,'b-');
% hold off
% hold on
% plot(Ti,cpi,'m-');
% hold off
% hold on
% plot(Ti,cnbi,'y-');  % (r,nQ,T) cost for Nbin
% hold off
% plot deterministic
sSTdetOptGraph_T(tmin,tmax,Kr,K0,mi,h,p,tstep);
