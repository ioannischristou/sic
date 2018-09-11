function z = sSTCIPpoissonGraph(ymin,ymax,T,L,lamda,h,p,p2)
if nargin < 8
    p2 = 0;
end
ylen = ymax-ymin+1;
y = 1:ylen;
G=1:ylen;
m=lamda*L;
for i=1:ylen
    y(i)= ymin+i-1;
    G(i) = sSTCIPpoisson(y(i),T,L,lamda,m,h,p,p2);
end
hold on
plot(y,G);
hold off
z = G;
end
