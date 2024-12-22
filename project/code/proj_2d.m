R = 1;                          
T = 0.2;             
k_const = 4; 
C_const = 2; 
f_func = @(r, t) 100 * exp(-10 * (r - 0.5).^2) * sin(2 * pi * t);

h = 0.025;                    
k = 0.1 * h^2;         
Nr = round(R / h) + 1;
Nt = round(T / k);
r = linspace(0, R, Nr)';

p = zeros(Nr, 1);     

p(1) = 0;   
p(end) = 0;

A = zeros(Nr, Nr);     
b = zeros(Nr, 1);      

for i = 2:Nr-1
    r_i = r(i);
    r_ip = (r(i) + r(i+1)) / 2;
    r_im = (r(i) + r(i-1)) / 2;

    A(i, i-1) = k_const * r_im / (r_i * h^2);
    A(i, i) = -k_const * (r_ip + r_im) / (r_i * h^2) - 1/k;
    A(i, i+1) = k_const * r_ip / (r_i * h^2);
end

A(1, 1) = -1/k - 2 * k_const / h^2;
A(1, 2) = 2 * k_const / h^2;

% Dirichlet
A(end, end) = 1;
A(end, end-1) = 0;
b(end) = 0;
%Neumann and Robin
%A(end, end) = alpha - beta / h;
%A(end, end-1) = beta / h;
%b(end) =g ; 

theta = linspace(0, 2*pi, 100); 
[R_grid, Theta_grid] = meshgrid(r, theta);
X = R_grid .* cos(Theta_grid);
Y = R_grid .* sin(Theta_grid);
P_plot = zeros(size(Theta_grid)); 

figure;
 shading interp;
colormap('hot');
colorbar;
caxis([0, 100]); 
xlabel('X');
ylabel('Y');
axis equal;

for n = 1:Nt
    b(2:Nr-1) = -p(2:Nr-1) / k - f_func(r(2:Nr-1), n*k);
    b(1) = -p(1) / k;        
    
    p = A \ b;
    
    P_plot = repmat(p', size(Theta_grid, 1), 1);
    pcolor(X, Y, P_plot);
    colorbar;
    shading interp;
    title(['Time = ', num2str(n*k)]);
    drawnow; 
end