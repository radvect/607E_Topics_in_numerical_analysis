format long g
t0 = 0;            
t_end = 1;           
y0 = 3;           
exact_solution = @(t) 12  ./ (exp(4 * t) + 3); 


f = @(t, y) y^2 - 4*y;


k_values = [0.1, 0.05, 0.01, 0.005];
errors = zeros(size(k_values)); 

for i = 1:length(k_values)
    k = k_values(i);
    t = t0:k:t_end;
    y = zeros(size(t));
    y(1) = y0;

    
    for n = 1:(length(t)-1)
        k1 = k * f(t(n), y(n));
        k2 = k * f(t(n) + k/2, y(n) + k1/2);
        k3 = k * f(t(n) + k/2, y(n) + k2/2);
        k4 = k * f(t(n) + k, y(n) + k3);
        y(n+1) = y(n) + (k1 + 2*k2 + 2*k3 + k4) / 6;
    end
    

    exact_y_end = exact_solution(t_end);
    errors(i) = abs(y(end) - exact_y_end);
end
p = polyfit(log(k_values), log(errors), 1);

figure;
loglog(k_values, errors, '-o', 'DisplayName', ['(Slope = ', num2str(p(1), '%.2f'), ')'], 'LineWidth', 1.5);
hold on
xlabel('Time Step k');
ylabel('Error at t = 1');
title('Convergence Study of RK4 Method');
grid on;
legend('show', 'Location', 'southeast');

hold off;

exportgraphics(gcf, "4b.png", 'Resolution', 300);
