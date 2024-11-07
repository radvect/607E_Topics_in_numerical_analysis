h_values = [0.1, 0.05, 0.01, 0.005, 0.0025, 0.00125]; 
errors = zeros(size(h_values)); 

for i = 1:length(h_values)
    h = h_values(i);

    x = (h:h:(1-h))';
    m = length(x);
    e = ones(m, 1);

    Dxx = spdiags([e (-2*e - h*h*e) e], [-1 0 1], m, m);
    Dxx = Dxx / h^2;


    f = 3.^x * (log(3)*log(3)-1)
    f(1) = f(1) - 1 / h^2;
    f(m) = f(m) - 3 / h^2;

    u_numerical = Dxx \ f;


    u_analytical = 3.^x;

    errors(i) = max(abs(u_numerical - u_analytical));
end


figure;
loglog(h_values, errors, '-o');
xlabel('Grid step h');
ylabel('Error');
title('Convergence Study for u_{xx} -u = 3^{x}(ln(x)^2-1)');
grid on;


p = polyfit(log(h_values), log(errors), 1);
hold on;
loglog(h_values, exp(polyval(p, log(h_values))), '--', 'DisplayName', sprintf('Slope: %.2f', p(1)));
legend('Error', sprintf('Slope: %.2f', p(1)), 'Location', 'best');
hold off;
exportgraphics(gcf, "5b.png", 'Resolution', 300);