
u = @(x) sin(x);
u_xxx_exact = @(x) -cos(x);

x0 = 1;
h_values = [0.1, 0.05, 0.01, 0.005];

error_5_point = zeros(size(h_values));
error_4_point = zeros(size(h_values));
error_6_point = zeros(size(h_values));

for i = 1:length(h_values)
    h = h_values(i);
    
    a_5 = (-1) / (2 * h^3);
    b_5 = 1 / h^3;
    c_5 = (-1) / h^3;
    d_5 = 1 / (2 * h^3);
    
  
    u_xxx_5_point = a_5 * u(x0 - 2*h) + b_5 * u(x0 - h) + c_5 * u(x0 + h) + d_5 * u(x0 + 2*h);
    error_5_point(i) = abs(u_xxx_5_point - u_xxx_exact(x0));
    

    a_4 = (-5) / (2 * h^3);
    b_4 = (18) / (2 * h^3);
    c_4 = (-24) / (2 * h^3);
    d_4 = (14) / (2 * h^3);
    e_4 = (-3) / (2 * h^3); 
    

    u_xxx_4_point = a_4 * u(x0) + b_4 * u(x0 + h) + c_4 * u(x0 + 2*h) + d_4 * u(x0 + 3*h) + e_4 * u(x0 + 4*h);
    error_4_point(i) = abs(u_xxx_4_point - u_xxx_exact(x0));
    

A = [
    1, 1, 1, 1, 1, 1;
    0, 1, 2, 3, 4, 5;
    0, 1/2, 2^2/2, 3^2/2, 4^2/2, 5^2/2;
    0, 1/6, 2^3/6, 3^3/6, 4^3/6, 5^3/6;
    0, 1/24, 2^4/24, 3^4/24, 4^4/24, 5^4/24;
    0, 1/120, 2^5/120, 3^5/120, 4^5/120, 5^5/120;
];
    b = [0; 0; 0; 1/(h^3); 0; 0];
    

    coeffs_6 = A \ b;
    a_6 = coeffs_6(1);
    b_6 = coeffs_6(2);
    c_6 = coeffs_6(3);
    d_6 = coeffs_6(4);
    e_6 = coeffs_6(5);
    f_6 = coeffs_6(6);

    u_xxx_6_point = a_6 * u(x0) + b_6 * u(x0 + h) + c_6 * u(x0 + 2*h) + d_6 * u(x0 + 3*h) + e_6 * u(x0 + 4*h) + f_6 * u(x0 + 5*h);
    error_6_point(i) = abs(u_xxx_6_point - u_xxx_exact(x0));
end

log_h = log(h_values);
slope_5_point = polyfit(log_h, log(error_5_point), 1);
slope_4_point = polyfit(log_h, log(error_4_point), 1);
slope_6_point = polyfit(log_h, log(error_6_point), 1);


figure;
loglog(h_values, error_4_point, '--s', 'DisplayName', ['4-Point (Slope = ', num2str(slope_4_point(1), '%.2f'), ')'], 'LineWidth', 1.5);
hold on;
loglog(h_values, error_5_point, '-o', 'DisplayName', ['5-Point (Slope = ', num2str(slope_5_point(1), '%.2f'), ')'], 'LineWidth', 1.5);
loglog(h_values, error_6_point, '-.d', 'DisplayName', ['6-Point (Slope = ', num2str(slope_6_point(1), '%.2f'), ')'], 'LineWidth', 1.5);
hold off;

xlabel('Grid step h');
ylabel('Error');
title('Convergence Study for 4,5,6 Points Schemes');
legend('show', 'Location', 'southeast');
grid on;
exportgraphics(gcf, "1d", 'Resolution', 300);