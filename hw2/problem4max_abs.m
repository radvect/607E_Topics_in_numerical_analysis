function demo_03_runge

figure(1); clf;
lw = 'linewidth'; ms = 'markersize';
lspace = linspace(-1.,1.,10000);

n = 8;

x_uni = linspace(-1,1,n+1)
x_cheb = cos(linspace(-pi, 0, n+1))
length(x_cheb)
%y_uni = 1./(1+x_uni.*x_uni);
%y_cheb = 1./(1+x_cheb.*x_cheb);

y_poly_uni = polynome(x_uni,lspace);
y_poly_cheb = polynome(x_cheb,lspace);


hold on
%original_line=plot(lspace,1./(1+lspace.^2),'-k', lw,2);
plot(x_uni, -0*ones(size(x_uni)), 'co', ms,12, lw,2);
%plot(x_uni, 1 ./ (1+x_uni.*x_uni), 'co', lw,2);
plot(x_cheb, -0*ones(size(x_cheb)), 'rx', ms,12, lw,2);
%plot(x_cheb, 1 ./ (1+x_cheb.*x_cheb), 'rx', lw,2);
h_uni_line = plot(lspace,y_poly_uni,'c', lw,2);
h_cheb_line = plot(lspace,y_poly_cheb,'r', lw,2);
%err = max(t1 - t');
%axis([-1.5 1.5])
xlim([-1.5, 1.5]);
title(sprintf('n = %2.0f', n))
xlabel('x')
ylabel('y')
yline(0, 'k--', 'LineWidth', 1); 
legend('pts uniform', 'pts chebyshev','pol uniform', 'pol chebyshev')
hold off

saveas(gcf, 'C:\\Users\\Pavel\\Documents\\MATLAB\\cheb.png'); 


function yi = polynome(x, xi)
  yi = prod(xi - x(:), 1); 
end


end
