L = 30;         
Tf = 17;          
h = 0.05;         
x = 0:h:L-h;        
N = length(x);      
nu = 1.1;%mb 0.8
k = nu * h;

u0 = exp(-20 * (x - 2).^2) + exp(-(x - 5).^2);  
u1 = exp(-20 * (x - k - 2).^2) + exp(-(x - k - 5).^2);

main_diag = zeros(N, 1);        
upper_diag = -nu * ones(N, 1);   
lower_diag = nu * ones(N, 1);  

A = spdiags([lower_diag, main_diag, upper_diag], [-1, 0, 1], N, N);
A(1, N) = nu;    
A(N, 1) = -nu;   

B = eye(N);

u_prev = u1';
u_prev_prev = u0';
target_time = 10;
for i = 1:round(Tf/k)
    u_next = A * u_prev + B * u_prev_prev;
    buffer = u_prev;
    u_prev = u_next;
    u_prev_prev = buffer;
    
    
    t = i * k;
    
    if abs(t - target_time) < k/2  
        figure;
        plot(x, u_next);
        axis([0 L -0.5 1.5]);
        title(sprintf('Solution at Time = %.2f', t));
        
   
        saveas(gcf, sprintf('3.png', t));
        break; 
    end
end

