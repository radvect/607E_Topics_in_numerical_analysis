function demo_03_runge
format long
lw = 'linewidth'; 
ms = 'markersize';
lspace = linspace(-1., 1., 10000);

max_uni = zeros(1, 29);  
max_cheb = zeros(1, 29); 

for n = 2:1:30
    x_uni = linspace(-1, 1, n);
    x_cheb = cos(linspace(-pi, 0, n));
    y_poly_uni = polynome(x_uni, lspace);
    y_poly_cheb = polynome(x_cheb, lspace);
    max_uni(n-1) = max(abs(y_poly_uni));
    max_cheb(n-1) = max(abs(y_poly_cheb));
end

n = 8;
x_uni = linspace(-1,1,n);
x_cheb = cos(linspace(-pi, 0, n));
y_poly_uni = polynome(x_uni,lspace);
y_poly_cheb = polynome(x_cheb,lspace);

%original_line=plot(lspace,1./(1+lspace.^2),'-k', lw,2);
%plot(x_uni, -0*ones(size(x_uni)), 'co', ms,12, lw,2);
%plot(x_uni, 1 ./ (1+x_uni.*x_uni), 'co', lw,2);
%plot(x_cheb, -0*ones(size(x_cheb)), 'rx', ms,12, lw,2);
semilogy(2:1:30, max_cheb, 'rx', ms, 12, lw, 2);
  
hold on;
xlabel('n');
ylabel('max(abs)');

semilogy(2:1:30, max_uni, 'co', ms, 12, lw, 2);

n_values = 2:1:30;

coeff_cheb = polyfit(2:1:30, log(max_cheb), 1); 
fit_cheb = polyval(coeff_cheb, 2:1:30);
semilogy(2:1:30, exp(fit_cheb), 'g-', lw, 2);  
pi_1 = (2 ./ n_values).^(n_values+1) .* factorial(n_values) / 4;
semilogy(2:1:30, (pi_1), 'b--', lw, 2); 

legend_text = sprintf('Experimental Bound cheb');
legend_text2 = sprintf('Theoretical Bound uniform');
legend( 'chebyshev pol','uniform pol', legend_text,legend_text2, 'Location', 'southwest');

saveas(gcf, 'C:\\Users\\Pavel\\Documents\\MATLAB\\bounds.png');
hold off;
saveas(gcf, 'C:\\Users\\Pavel\\Documents\\MATLAB\\bounds.png');

function yi = polynome(x, xi)
    yi = prod(xi - x(:), 1); 
end

end
