format long g
t0 = 0;
t_end = 1;
y0 = 3;
exact_solution = @(t) 12 ./ (exp(4 * t) + 3);

A = @(y) y^2;
B = @(y) -4*y;

k_start = 0.1;
N = 20;
k_values = k_start ./ (2.^(0:N-1));
errors = zeros(size(k_values));

for i = 1:length(k_values)
    k = k_values(i);
    t = t0:k:t_end;
    y = zeros(size(t));
    y(1) = y0;

    for n = 1:(length(t)-1)

        k1_A = A(y(n));
        k2_A = A(y(n) + 0.5 * k * k1_A);
        k3_A = A(y(n) + 0.5 * k * k2_A);
        k4_A = A(y(n) + k * k3_A);
        yn = y(n) + (k / 6) * (k1_A + 2*k2_A + 2*k3_A + k4_A);


        k1_B = B(yn);
        k2_B = B(yn + 0.5 * k * k1_B);
        k3_B = B(yn + 0.5 * k * k2_B);
        k4_B = B(yn + k * k3_B);
        y(n+1) = yn + (k / 6) * (k1_B + 2*k2_B + 2*k3_B + k4_B);
    end

    exact_y_end = exact_solution(t_end);
    errors(i) = abs(y(end) - exact_y_end);
end

p = polyfit(log(k_values), log(errors), 1);

figure;
loglog(k_values, errors, '-o', 'DisplayName', ['Operator Splitting (RK4 Slope = ', num2str(p(1), '%.2f'), ')'], 'LineWidth', 1.5);
hold on
xlabel('Time Step k');
ylabel('Error at t = 1');
title('Convergence Study of Operator Splitting with RK4');
grid on;
legend('show', 'Location', 'southeast');

% Save the figure
exportgraphics(gcf, "4c.png", 'Resolution', 300);

hold off;