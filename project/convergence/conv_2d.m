function main()
    R = 1;
    Tmax = 2;
    alpha = 1;
    f_func = @(r, t) 100 * exp(-10 * (r - 0.5).^2) * sin(2 * pi * t); 

    Nr_values = [20, 40, 80, 160, 320];
    Nt = 200;

    ref_Nr = 3000;
    ref_u = compute_reference_solution(R, Tmax, alpha, f_func, ref_Nr, Nt);

    h_values = R ./ (Nr_values - 1);
    errors = compute_errors(Nr_values, Nt, R, Tmax, alpha, f_func, ref_u);

    fit_and_plot(h_values, errors);
end

function ref_u = compute_reference_solution(R, Tmax, alpha, f_func, ref_Nr, Nt)
    ref_dr = R / (ref_Nr - 1);
    ref_r = linspace(0, R, ref_Nr)';
    ref_u = zeros(ref_Nr, 1);
    ref_dt = Tmax / Nt;

    A_ref = assemble_matrix(ref_Nr, ref_r, alpha, ref_dr, ref_dt);
    
    for n = 1:Nt
        b_ref = assemble_rhs(ref_Nr, ref_u, ref_r, f_func, ref_dt, n);
        ref_u = A_ref \ b_ref;
    end
end

function errors = compute_errors(Nr_values, Nt, R, Tmax, alpha, f_func, ref_u)
    ref_Nr = length(ref_u);
    ref_r = linspace(0, R, ref_Nr)';
    ref_dr = R / (ref_Nr - 1);

    errors = zeros(length(Nr_values), 1);

    for iNr = 1:length(Nr_values)
        Nr = Nr_values(iNr);
        dr = R / (Nr - 1);
        r = linspace(0, R, Nr)';
        u = zeros(Nr, 1);
        dt = Tmax / Nt;

        A = assemble_matrix(Nr, r, alpha, dr, dt);

        for n = 1:Nt
            b = assemble_rhs(Nr, u, r, f_func, dt, n);
            u = A \ b;
        end

        u_interp = interp1(r, u, ref_r, 'linear', 'extrap');
        errors(iNr) = sqrt(sum((u_interp - ref_u).^2) * ref_dr);
    end
end

function fit_and_plot(h_values, errors)
    log_h = log(h_values);
    log_errors = log(errors);
    p = polyfit(log_h, log_errors, 1);
    log_errors_fit = polyval(p, log_h);

    figure;
    loglog(h_values, errors, '-o', 'DisplayName', 'Numerical Error');
    hold on;
    loglog(h_values, exp(log_errors_fit), '--', 'DisplayName', sprintf('Fit: slope = %.2f', p(1)));
    xlabel('Step size (h)');
    ylabel('L2 Error');
    title('Convergence Study');
    legend('Location', 'best');
    grid on;
end

function A = assemble_matrix(Nr, r, alpha, dr, dt)
    A = zeros(Nr, Nr);

    for i = 2:Nr-1
        r_i = r(i);
        r_ip = (r(i) + r(i+1)) / 2;
        r_im = (r(i) + r(i-1)) / 2;

        A(i, i-1) = alpha * r_im / (r_i * dr^2);
        A(i, i) = -alpha * (r_ip + r_im) / (r_i * dr^2) - 1/dt;
        A(i, i+1) = alpha * r_ip / (r_i * dr^2);
    end

    A(1, 1) = -1/dt - 2 * alpha / dr^2;
    A(1, 2) = 2 * alpha / dr^2;
    A(end, end) = 1;
    A(end, end-1) = 0;
end

function b = assemble_rhs(Nr, u, r, f_func, dt, time_step)
    b = zeros(Nr, 1);
    b(2:Nr-1) = -u(2:Nr-1) / dt - f_func(r(2:Nr-1), time_step * dt);
    b(1) = -u(1) / dt;
    b(end) = 0;
end