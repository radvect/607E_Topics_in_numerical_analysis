% Demo01
% Approximately solve the heat eqn u_t = u_xx
% on 0 <= x < 1, with IC u0 = sin(2*pi*x)
% with periodic BCs, on 0 <= t <= Tf.

%% choose a spatial discretization
h = 1/50
x = 0:h:1; x = x'
N = length(x)

u0 = sin(2*pi*x)

%figure(1); clf;
%plot(x, u0, 'rx')
%hold on
%xlabel('x'); ylabel('u')

%% Build a matrix that approximates the second deriv operator
e = ones(size(x));
L = spdiags([e -2*e e], [-1 0 1], N, N);
%L(1,N) = 1;
%L(N,1) = 1;
L(1,1) = 1;
L(1,2)=0;
L(N,N-1)=0;
L(N,N) = 1;
%% to look at the matrix directly use
%full(L)


L = (1/h^2)*L;
full(L*h^2)

%% testing the 2nd derivative
% uxx_exact = -4*pi^2*sin(2*pi*x);
%uxx = L*u0;
% plot(x, uxx, 'b-o')
% plot(x, uxx_exact, 'gx')
% max(uxx - uxx_exact)

Tf = 0.02;
k = 0.25*h*h;
numsteps = ceil(Tf / k);
k = Tf / numsteps;

u = u0;
handle = plot(x, u0, 'k-o');
legend('u0(x)', 'u(x,t)')
for n=1:numsteps
  unew = u + k*(L*u);
  set(handle, 'ydata', unew)
  u = unew;
  title(['t = ' num2str(n*k)])
  drawnow
end
pause
