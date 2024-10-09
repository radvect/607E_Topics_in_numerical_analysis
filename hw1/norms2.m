Tf = 0.25;
%h_values = [0.125, 0.0625, 0.03125, 0.015625]; % Массив значений h
h_values = 2.^(-(2:10));
L2_final_results = zeros(1, length(h_values)); % Массив для хранения L2 для последней итерации каждого h

for idx = 1:length(h_values)
    %idx
    h = h_values(idx);
    h
    dt = 0.25 * h * h;
    dt
    x = 0:h:(1-h);
    x = x';
    N = length(x);
    e = ones(size(x));

    L = spdiags([e -2*e e], [-1 0 1], N, N);
    L(1, N) = 1;
    L(N, 1) = 1;
    numsteps = ceil(Tf / dt);
    L = (1 / h^2) * L;
    u0 = sin(2 * pi * x); % Начальное условие
    u = u0;

    for n = 0:numsteps
        curr_t = n * dt;
        X_exact_t = exp(-4 * pi*pi * curr_t) * sin(2 * pi * x); % Точное решение
        L2_value = max_norm = norm(X_exact_t - u, inf);; % Вычисление нормы L2

        % На последней итерации сохранить L2_value
        if n == numsteps
	    %norm(X_exact_t)
	    %L
	    L2_final_results(idx) = L2_value;
        end

        unew = u + dt * (L * u);
        u = unew;
    end
end


h_values
L2_final_results


% Assuming L2_final_results and h_values are already calculated from your previous code
log_h = log(h_values); % Take the logarithm of h values
log_L2 = log(L2_final_results); % Take the logarithm of L2 norm values

% Use polyfit to find the best-fit line (linear fit in log-log space)
p = polyfit(log_h, log_L2, 1); % Linear fit in log-log space

% Generate a set of points for the best-fit line
fit_L2 = polyval(p, log_h);

% Plot the log-log plot with data points
figure;
loglog(h_values, L2_final_results, 'bo-', 'LineWidth', 2); % Plot the actual data
hold on;
loglog(h_values, exp(fit_L2), 'r--', 'LineWidth', 2); % Plot the best-fit line
hold off;

% Add labels and title
xlabel('h');
ylabel('Error');
title('Log-Log Plot of Error vs h');
legend('Error', 'Best Fit Line', 'Location', 'SouthEast');
grid on;
slope = p(1);

text_pos_x = min(h_values) * 1.5; % Position for the text (adjust as needed)
text_pos_y = max(L2_final_results) * 0.5; % Position for the text (adjust as needed)

text(text_pos_x, text_pos_y, sprintf('Slope: %.2f', slope), 'FontSize', 12, 'Color', 'black', 'BackgroundColor', 'white');

% Display the slope of the line
disp(['Best-fit line slope: ', num2str(p(1))]);
pause
