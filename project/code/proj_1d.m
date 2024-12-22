h = 0.025;                    
k = 0.5 * h^2;                 
x = (0 : h : 1)';            
N = length(x);


p0 = @(x) cos(2*pi*x);  
p = p0(x);                     


f = @(x, t) cos(2*pi*x)*exp(-t);    

alpha_left = 0; beta_left = 1; gamma_left = 0;  
alpha_right = 0; beta_right = 1; gamma_right = 0; 


a = -2 / h^2;
b = 1 / h^2;
main_diag = a * ones(N, 1);
off_diag = b * ones(N-1, 1);
L = diag(main_diag) + diag(off_diag, 1) + diag(off_diag, -1);

%Dirichlet
%Nothing should be changed.
alpha = 0
beta = 1
%Robin and Neumann 
L(1,1) = -2/h/h - 2 *alpha/beta/h
L(1,2) = 2/h/h
L(end,end-1) = 2/h/h
L(end,end) = -2/h/h - 2 *alpha/beta/h


Tf = 0.5;                
numsteps = ceil(Tf / k);


figure(1); clf;
plt = plot(x, p, 'k.-', 'linewidth', 2);
axis([0 1 -1.1 1.1])
grid on
xlabel('x'); ylabel('p(t)');

k_const = 5 / (4 * pi^2);
C_const = 4 ;


t = 0;                     
for n = 1:numsteps
    f_values = f(x, t);
    p_new = p + k/C_const* (k_const *L * p - f_values);
    p = p_new;


    t = t + k;
    set(plt, 'YData', p)
    set(title(sprintf('t = %.3f', t)), 'Interpreter', 'none');
    drawnow
end