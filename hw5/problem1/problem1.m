h_values = 0.1 ./ (2.^(0:10-1)); 
errors = zeros(size(h_values)); 

for i = 1:length(h_values)
    h = h_values(i);

    x = (0:h:1)';
    m = length(x);
    e = ones(m, 1);
    
    Dxx = spdiags([e (-2*e - h*h*e) e], [-1 0 1], m, m);

    Dxx(1,1) = -(2+h*h);
    Dxx(m,m) = -2-h*h;
    Dxx(1,2) = 2
    Dxx(m, m-1)=0
    %disp(Dxx)
    
    Dxx = Dxx / h^2;

    f = (-4*pi*pi-1)*cos(2*pi*x)%3.^x * (log(3)*log(3)-1);
    f(1) = f(1)
    f(m) = f(m) - 2/h^2
    u_numerical = Dxx \ f;


    u_analytical = cos(2*pi*x)%3.^x;
    
    errors(i) = max(abs(u_numerical - u_analytical));
end


figure;
loglog(h_values, errors, '-o');
xlabel('Grid step h');
ylabel('Error');
title('Convergence Study for u_{xx} -u = -cos(2\pi x)(4\pi^2+1)');
grid on;


p = polyfit(log(h_values), log(errors), 1);
hold on;
loglog(h_values, exp(polyval(p, log(h_values))), '--', 'DisplayName', sprintf('Slope: %.2f', p(1)));
legend('Error', sprintf('Slope: %.2f', p(1)), 'Location', 'best');
hold off;
exportgraphics(gcf, "convergence_study.png", 'Resolution', 300);
