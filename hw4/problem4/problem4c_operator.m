


t0 = 0;               
t_end = 2;          
y0 = 3;               
k = 0.01;           
exact_solution = @(t) 12  ./ (exp(4 * t) + 3); 

A = @(y) y^2;     
B = @(y) -4*y; 


t = t0:k:t_end;
y = zeros(size(t));
y(1) = y0;

for n = 1:(length(t)-1)

    yn = y(n) + k * A(y(n));      
    y(n+1) = yn + k * B(yn);   
end


figure;
plot(t, y, 'b-', 'DisplayName', 'Numerical Solution (RK4 with Operator Splitting)');
hold on;
plot(t, exact_solution(t), 'r--', 'DisplayName', 'Exact Solution');
xlabel('Time t');
ylabel('y');
title('y''(t) = y^2 - 4y; RK4 with Operator Splitting');
legend('Location', 'best');
grid on;

error_at_t1 = abs(y(end) - exact_solution(t_end));
disp(['t = 1: ', num2str(error_at_t1)]);
exportgraphics(gcf, "4c_sol.png", 'Resolution', 300);