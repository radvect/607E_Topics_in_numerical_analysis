h_values = 0.1 ./ (2.^(0:10-1)); 
errors = zeros(size(h_values)); 

for i = 1:length(h_values)
    h = h_values(i);

    x = (0:h:(1))';
    m = length(x);
    e = ones(m, 1);

    Dxx = spdiags([e (-2*e - h*h*e) e], [-1 0 1], m, m)
    Dxx(1,1) = -(2+h*h);
    Dxx(m,m) = -(2+h*h);
    Dxx(1,2) =0;
    Dxx(m, m-1)=0;
    Dxx = Dxx / h^2;


    f = -sin(pi * x) * (pi^2 + 1);

    u_numerical = Dxx \ f;


    u_analytical = sin(pi * x);

    errors(i) = max(abs(u_numerical - u_analytical));
end


figure;
loglog(h_values, errors, '-o');
xlabel('Grid step h');
ylabel('Error');
title('Convergence Study for u_{xx} -u = -sin(\pi x)(\pi^2+1)');
grid on;

% Optional: Fit a line to check convergence rate
p = polyfit(log(h_values), log(errors), 1);
hold on;
loglog(h_values, exp(polyval(p, log(h_values))), '--', 'DisplayName', sprintf('Slope: %.2f', p(1)));
legend('Error', sprintf('Slope: %.2f', p(1)), 'Location', 'best');
hold off;
exportgraphics(gcf, "6d1.png", 'Resolution', 300);
