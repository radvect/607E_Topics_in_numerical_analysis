

e = 1/100;
fun = @(x,y) exp(-((x - 1.5).^2+(y - 0.6).^2)/(1/100));

xmin = 0;
xmax =2;
ymin = 0;
ymax = 1;

q = integral2(fun,xmin,xmax,ymin,ymax)


