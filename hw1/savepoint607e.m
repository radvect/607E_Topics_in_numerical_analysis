% Demo01
% Approximately solve the heat eqn u_t = u_xx
% on 0 <= x < 1, with IC u0 = sin(2*pi*x)
% with periodic BCs, on 0 <= t <= Tf.

%% choose a spatial discretization
h = 1/50
x = 0:h:(1-h); x = x'
N = length(x)

u0 = sin(2*pi*x)

%figure(1); clf;
%plot(x, u0, 'rx')
%hold on
%xlabel('x'); ylabel('u')

%% Build a matrix that approximates the second deriv operator
e = ones(size(x));
L = spdiags([e -2*e e], [-1 0 1], N, N);
L(1,N) = 1;
L(N,1) = 1;
%% to look at the matrix directly use
% full(L)
L = (1/h^2)*L;

%% testing the 2nd derivative
% uxx_exact = -4*pi^2*sin(2*pi*x);
% uxx = L*u0;
% plot(x, uxx, 'b-o')
% plot(x, uxx_exact, 'gx')
% max(uxx - uxx_exact)


Tf = 0.25;
dt = 0.25*h*h;
numsteps = ceil(Tf / dt);
%k = Tf / numsteps;
t = 0:dt:Tf; t = t'
[T, X] = meshgrid(t, x);

X_exact = exp(-2*pi*T) .* sin(2*pi*X);
L2 = []; 
u = u0;
%handle = plot(x, u0, 'k-o');
%legend('u0(x)', 'u(x,t)')

for n=0:numsteps
  curr_t = n*dt
  X_exact_t = exp(-2*pi*curr_t)*sin(2*pi*x)  
  L2_value = norm(X_exact_t- u)
  unew = u + dt*(L*u);
  L2(end+1)=L2_value
  %set(handle, 'ydata', unew)
u = unew;

  %title(['t = ' num2str(n*dt)])
  %drawnow
end

plot(t, L2)
xlabel('X-axis'); % Название оси x
ylabel('Y-axis'); % Название оси y
title('Plot of X vs Y'); % Заголовок графика
grid on; % Включить сетку
pause
drawnow
