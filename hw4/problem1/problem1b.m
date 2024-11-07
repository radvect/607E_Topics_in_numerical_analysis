

u = @(x) sin(x);        
u_x_exact = @(x) cos(x); 

x0 = 1;


h_values = [0.1, 0.05, 0.01, 0.005];


error_forward = zeros(size(h_values));
error_backward = zeros(size(h_values));
error_centered = zeros(size(h_values));
error_one_sided = zeros(size(h_values));

for i = 1:length(h_values)
    h = h_values(i);
    
    u_x_forward = (u(x0 + h) - u(x0)) / h;
    error_forward(i) = abs(u_x_forward - u_x_exact(x0));
    
    
    u_x_backward = (u(x0) - u(x0 - h)) / h;
    error_backward(i) = abs(u_x_backward - u_x_exact(x0));

    u_x_centered = (u(x0 + h) - u(x0 - h)) / (2 * h);
    error_centered(i) = abs(u_x_centered - u_x_exact(x0));

    u_x_one_sided = (-3/2 * u(x0) + 2 * u(x0 + h) - 1/2 * u(x0 + 2*h)) / h;
    error_one_sided(i) = abs(u_x_one_sided - u_x_exact(x0));
end


log_h = log(h_values);
slope_forward = polyfit(log_h, log(error_forward), 1);
slope_backward = polyfit(log_h, log(error_backward), 1);
slope_centered = polyfit(log_h, log(error_centered), 1);
slope_one_sided = polyfit(log_h, log(error_one_sided), 1);

figure;
loglog(h_values, error_forward, '-o', 'DisplayName', ['Forward (Slope = ', num2str(slope_forward(1), '%.2f'), ')'], 'LineWidth', 1.5);
hold on;
loglog(h_values * 1.05, error_backward, '--s', 'DisplayName', ['Backward (Slope = ', num2str(slope_backward(1), '%.2f'), ')'], 'LineWidth', 1.5); % slight shift
loglog(h_values, error_centered, '-.^', 'DisplayName', ['Centered (Slope = ', num2str(slope_centered(1), '%.2f'), ')'], 'LineWidth', 1.5);
loglog(h_values, error_one_sided, '-.d', 'DisplayName', ['One-Sided (Slope = ', num2str(slope_one_sided(1), '%.2f'), ')'], 'LineWidth', 1.5);
hold off;

xlabel('h');
ylabel('Error');
title('Convergence Study Test for u=sin(x), x=1');
legend('show', 'Location', 'southeast');
grid on;
exportgraphics(gcf, "1b.png", 'Resolution', 300);