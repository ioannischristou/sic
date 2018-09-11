function [sopti Sopti Topti copti timesi] = sSTCpoissonOpt_MultiKr(Tmax, Kri, K0, L, lamda, h, p, p2, epst, tnot, estop)
% Kri(1) must be largest element in vector
cases = numel(Kri);
sopti = 1:cases;
Sopti = 1:cases;
Topti = 1:cases;
copti = 1:cases;
timesi = 1:cases;
disp(['calling sSTCpoissonOpt() for Kr=' num2str(Kri(1)) ]);
tic;
[s1 S1 T1 c1 s1i S1i flen T1i] = sSTCpoissonOpt(Tmax, Kri(1), K0, L, lamda, h, p, p2, epst, tnot, estop);
timesi(1)=toc;
sopti(1)=s1;
Sopti(1)=S1;
Topti(1)=T1;
copti(1)=c1;
for i=2:cases
   disp(['Computing Opt. for Kr=' num2str(Kri(i))]);
   timesi(i)=timesi(1);
   cbi=10^30;
   for j=1:flen
       c = sSTCpoisson(s1i(j), S1i(j), T1i(j), Kri(i), K0, L, lamda, h, p, p2);
       if c < cbi
           sopti(i)=s1i(j);
           Sopti(i)=S1i(j);
           Topti(i)=T1i(j);
           copti(i)=c;
           cbi = c;
       end
   end
end

end

