Tf = 0.25;
h_values = [0.005]; % Array of h values
L2_final_results = zeros(1, length(h_values)); % Array to store L2 norms for each h

for idx = 1:length(h_values)
    h = h_values(idx);
    dt = 0.25 * h * h;
    x = 0:h:(1 - h);
    x = x';
    N = length(x);
    e = ones(N, 1);

    % Create matrix L with periodic boundary conditions
    L = spdiags([e -2*e e], [-1 0 1], N, N);
    L(1, N) = 1;
    L(N, 1) = 1;
    L = (1 / h^2) * L;

    % Initial condition
    u0 = sin(pi * x);
    u = u0;

    numsteps = ceil(Tf / dt);

    % Initialize arrays to store solutions at specific times
    numerical_t001 = [];
    numerical_t01 = [];
    numerical_t02 = [];
    threeterm_t001 = [];
    threeterm_t01 = [];
    threeterm_t02 = [];
    twoterm_t001 = [];
    twoterm_t01 = [];
    twoterm_t02 = [];

    % Main time-stepping loop
    for n = 0:numsteps
        curr_t = n * dt;

        % Three-term analytical solution
        X_exact_t = (2/pi) + (4)/(-3*pi) * cos(2*pi*x) * exp(-4*pi^2*curr_t) ...
                    + (4)/(-15*pi) * cos(4*pi*x) * exp(-16*pi^2*curr_t);

        % Numerical solution at the current time step
        if n > 0
            u = u + dt * (L * u);
        end

        % Compute L-infinity norm of the error
        L2_value = norm(X_exact_t - u, inf);

        % Save results at specified time instances
        if abs(curr_t - 0.01) < dt/2
            threeterm_t001 = X_exact_t;
            twoterm_t001 = (2/pi) + (4)/(-3*pi) * cos(2*pi*x) * exp(-4*pi^2*curr_t);
            numerical_t001 = u;
            L2001 = L2_value;
        end
        if abs(curr_t - 0.1) < dt/2
            threeterm_t01 = X_exact_t;
            twoterm_t01 = (2/pi) + (4)/(-3*pi) * cos(2*pi*x) * exp(-4*pi^2*curr_t);
            numerical_t01 = u;
            L201 = L2_value;
        end
        if abs(curr_t - 0.2) < dt/2
            threeterm_t02 = X_exact_t;
            twoterm_t02 = (2/pi) + (4)/(-3*pi) * cos(2*pi*x) * exp(-4*pi^2*curr_t);
            numerical_t02 = u;
            L202 = L2_value;
        end
    end

    % Display norms of errors
    disp(['Error norm at t = 0.01: ', num2str(L2001)]);
    disp(['Error norm at t = 0.1: ', num2str(L201)]);
    disp(['Error norm at t = 0.2: ', num2str(L202)]);
end

% Plotting the results
figure;

% Plots for t = 0.01
subplot(2, 3, 1);
plot(x, numerical_t001, 'b-', 'LineWidth', 1.5); hold on;
plot(x, threeterm_t001, 'r--', 'LineWidth', 1.5);
title('t = 0.01: Numerical and Analytical (3-term)');
xlabel('x'); ylabel('u(x, t)');
legend('Numerical', 'Analytical (3-term)');
grid on;

subplot(2, 3, 4);
plot(x, abs(numerical_t001 - threeterm_t001), 'k-', 'LineWidth', 1.5); hold on;
plot(x, abs(numerical_t001 - twoterm_t001), 'g--', 'LineWidth', 1.5);
title('t = 0.01: Error');
xlabel('x'); ylabel('Error');
legend('Numerical - Analytical (3-term)', 'Numerical - Analytical (2-term)');
grid on;

% Plots for t = 0.1
subplot(2, 3, 2);
plot(x, numerical_t01, 'b-', 'LineWidth', 1.5); hold on;
plot(x, threeterm_t01, 'r--', 'LineWidth', 1.5);
title('t = 0.1: Numerical and Analytical (3-term)');
xlabel('x'); ylabel('u(x, t)');
legend('Numerical', 'Analytical (3-term)');
grid on;

subplot(2, 3, 5);
plot(x, numerical_t01 - threeterm_t01, 'k-', 'LineWidth', 1.5); hold on;
plot(x, numerical_t01 - twoterm_t01, 'g--', 'LineWidth', 1.5);
title('t = 0.1: Error');
xlabel('x'); ylabel('Error');
legend('Numerical - Analytical (3-term)', 'Numerical - Analytical (2-term)');
grid on;

% Plots for t = 0.2
subplot(2, 3, 3);
plot(x, numerical_t02, 'b-', 'LineWidth', 1.5); hold on;
plot(x, threeterm_t02, 'r--', 'LineWidth', 1.5);
title('t = 0.2: Numerical and Analytical (3-term)');
xlabel('x'); ylabel('u(x, t)');
legend('Numerical', 'Analytical (3-term)');
grid on;

subplot(2, 3, 6);
plot(x, numerical_t02 - threeterm_t02, 'k-', 'LineWidth', 1.5); hold on;
plot(x, numerical_t02 - twoterm_t02, 'g--', 'LineWidth', 1.5);
title('t = 0.2: Error');
xlabel('x'); ylabel('Error');
legend('Numerical - Analytical (3-term)', 'Numerical - Analytical (2-term)');
grid on;

% Overall title for all plots
%sgtitle('Comparison of Numerical and Analytical Solutions and Errors');
pause
