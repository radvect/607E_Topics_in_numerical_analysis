h = 0.01;

x = (h:h:(1-h))';
m = length(x);
e = ones(m, 1);

Dxx = spdiags([e (-2*e - h*h*e) e], [-1 0 1], m, m);
Dxx = Dxx / h^2;

f = -sin(pi * x) * (pi^2 + 1);

u_numerical = Dxx \ f;

u_analytical = sin(pi * x);

error = max(abs(u_numerical - u_analytical));

figure;
plot(x, u_numerical, '-o', 'DisplayName', 'Numerical Solution');
hold on;
plot(x, u_analytical, '-', 'DisplayName', 'Analytical Solution');
xlabel('x');
ylabel('u');
title(sprintf('Solution for u_{xx} - u = -sin(\\pi x)(\\pi^2+1) with h = %.4f', h));
legend;
grid on;
hold off;
exportgraphics(gcf, "5a_add.png", 'Resolution', 300);