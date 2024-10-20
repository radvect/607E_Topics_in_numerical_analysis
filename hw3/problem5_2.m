function [u11_values, u12_values] = solve_FE(U0, num_iter)
    k = 17.0652165601579625588917206249 / 29000;
    mu = 0.012277471;
    mu_bar = 1 - mu;
    u11_values = zeros(num_iter, 1);
    u12_values = zeros(num_iter, 1);
    U_n = U0;

    for i = 1:num_iter
        u11_values(i) = U_n(1);
        u12_values(i) = U_n(2);
        D1 = ((U_n(1) + mu)^2 + U_n(2)^2)^(3/2);
        D2 = ((U_n(1) - mu_bar)^2 + U_n(2)^2)^(3/2);

        A = [1, 0, k, 0;
             0, 1, 0, k;
             k * (1 - mu_bar / D1 - mu / D2), 0, 1, 2 * k;
             0, k * (1 - mu_bar / D1 - mu / D2), -2 * k, 1];
        B = [0; 0; k * mu * mu_bar * (-1 / D1 + 1 / D2); 0];
        U_next = A * U_n + B;
        U_n = U_next;
    end
end
function [t, U] = solve_ode45(U0, num_iter)
    k = 17.0652165601579625588917206249 / 29000;
    options = odeset('AbsTol', 1e-11, 'RelTol', 1e-8);
    t_span = [0, num_iter * k];

    function dUdt = three_body_ode(~, U)
        mu = 0.012277471;
        mu_bar = 1 - mu;
        u1 = U(1); u2 = U(2); v1 = U(3); v2 = U(4);
        D1 = ((u1 + mu)^2 + u2^2)^(3/2);
        D2 = ((u1 - mu_bar)^2 + u2^2)^(3/2);
        du1dt = v1;
        du2dt = v2;
        dv1dt = 2 * v2 + u1 - mu_bar * (u1 + mu) / D1 - mu * (u1 - mu_bar) / D2;
        dv2dt = -2 * v1 + u2 - mu_bar * u2 / D1 - mu * u2 / D2;
        dUdt = [du1dt; du2dt; dv1dt; dv2dt];
    end

    [t, U] = ode45(@three_body_ode, t_span, U0, options);
end

U0 = [0.994; 0.0; 0.0; -2.00158510637908252240537862224];
num_iter = 29000;
[t_ode, U_ode] = solve_ode45(U0, num_iter);
y1_ode45 = U_ode(:, 1);
y2_ode45 = U_ode(:, 2);
[y1_FE, y2_FE] = solve_FE(U0, num_iter);

figure;
plot(y1_FE, y2_FE, 'b', 'LineWidth', 1.5);
hold on;
plot(y1_ode45, y2_ode45, 'r', 'LineWidth', 1.5);
xlabel('y1'); ylabel('y2');
title('Trajectory of the spacecraft');
legend('FE', 'ode45');
grid on;
plot(-0.012277471, 0, 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
plot(1 - 0.012277471, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
legend('FE', 'ode45', 'Earth', 'Moon');
hold off;
saveas(gcf, 'C:\\Users\\Pavel\\Documents\\MATLAB\\hw3\\pic\\spacecraft_fe_ode45.png');

