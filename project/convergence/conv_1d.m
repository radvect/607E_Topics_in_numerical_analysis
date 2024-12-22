h_values = 0.1 ./ (2.^(0:3));  
errors = zeros(size(h_values));

Tf = 0.5;
k_factor = 1;

for idx = 1:length(h_values)
    h = h_values(idx);
    k = 0.01* k_factor * h^2
    
    %/2/(k_const/C);
    x = (h:h:1)';
    N = length(x);
    u0 = @(x) cos(2*pi*x);
    u = u0(x);
    f = @(x, t) -cos(2*pi*x)*exp(-t);

    a = -2 / h^2;
    b = 1 / h^2;
    main_diag = a * ones(N, 1);
    off_diag = 1 / h^2 * ones(N-1, 1);
    L = diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1);

    alpha = 0
    beta = 1

    L(1,1) = -2/h/h - 2 *alpha/beta/h
    L(1,2) = 2/h/h
    L(end,end-1) = 2/h/h
    L(end,end) = -2/h/h - 2 *alpha/beta/h


    k_const = 1 / (2 * pi^2);
    t = 0;
    C = 1; 
    numsteps = ceil(Tf / k);
    for n = 1:numsteps
        f_values = f(x, t);
        u_new = u + k * (k_const/C * L * u - f_values);
        u = u_new;
        t = t + k;
    end

    u_exact = cos(2*pi*x) * exp(-Tf);

    errors(idx) = max(abs(u - u_exact));
end

log_h = log(h_values);
log_errors = log(errors);
p = polyfit(log_h, log_errors, 1); 
log_errors_fit = polyval(p, log_h);
errors_fit = exp(log_errors_fit);


figure;
loglog(h_values, errors, '-o', 'LineWidth', 2);  
hold on;
loglog(h_values, errors_fit, '--', 'LineWidth', 2);  
xlabel('Step size (h)');
ylabel('L2 Error');
title('Convergence Study');
legend('Location', 'best');
legend('Numerical Error', sprintf('Fit: slope = %.2f', p(1)), 'Location', 'best');
grid on;

