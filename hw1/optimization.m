pkg load optim;
format long
r =1;

function draw_circle(x, y, r, color)
    theta = linspace(0, 2*pi, 100);
    x_circle = r * cos(theta) + x;
    y_circle = r * sin(theta) + y;
    plot(x_circle, y_circle, 'Color', color, 'LineWidth', 2);
    axis equal;
end

function cost = objective_function(vars)

    R = vars(end);
    x = vars(1:2:end-1);
    y = vars(2:2:end-1);
    N = length(x);

    cost = R*R
end


function [c, ceq] = constraints(vars)
    N = (length(vars) - 1) / 2;
    x = vars(1:2:end-1);
    y = vars(2:2:end-1);
    R = vars(end);
    r=1;

    c_inside = sqrt(x.^2 + y.^2) + r - R;


    c_overlap = [];
    for i = 1:N-1
        for j = i+1:N
            distance = sqrt((x(i) - x(j))^2 + (y(i) - y(j))^2);
            c_overlap(end+1) = 2*r - distance;
        end
    end

    c = [c_inside; c_overlap'];
    ceq = [];
end


N =11 ;


theta = linspace(0, 2*pi, N+1);
theta(end) = [];
x0 = zeros(2*N + 1, 1);
M = N
r_initial =  linspace(0, 10, M);;

for j = 1:M
  for i = 1:N
      x0(2*i -1) = r_initial(end) * cos(theta(i))%+5*rand(1);
      x0(2*i) = r_initial(end) * sin(theta(i)) %+5*rand(1);
end
end
x0(end) = r_initial(M) + 5*r;
figure;
hold on;
colors = hsv(N);

for i = 1:N
    draw_circle(x0(2*i -1), x0(2*i), r, colors(i,:));
    plot(x0(2*i -1), x0(2*i), 'o', 'Color', colors(i,:), 'MarkerSize', 10, 'MarkerFaceColor', colors(i,:));
end

draw_circle(0, 0, x0(end), 'g');  % Big circle
plot(0, 0, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
title('Optimized Circle Packing');
axis equal;
hold off;

pause;

global r;
[c, ceq] = constraints(x0);
disp('Constraint values for x0 (should be <= 0):');
disp(c);


if all(c <= 0)
    disp('Initial guess is feasible.');
else
    disp('Initial guess violates constraints.');
end

options = optimset('Algorithm', 'interior-point', 'Display', 'iter', "MaxIter",11);

[x_opt, fval] = fmincon(@objective_function, x0, [], [], [], [], [], [], @constraints, options);


disp('Optimized vector:');
disp(x_opt);
disp('Minimum radius R:');
disp(sqrt(fval));

figure;
hold on;
colors = hsv(N);

for i = 1:N
    draw_circle(x_opt(2*i -1), x_opt(2*i), r, colors(i,:));
    plot(x_opt(2*i -1), x_opt(2*i), 'o', 'Color', colors(i,:), 'MarkerSize', 10, 'MarkerFaceColor', colors(i,:));
end
draw_circle(0, 0, x_opt(end), 'g');
plot(0, 0, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
title('Optimized Circle Packing');
% Выносим значение минимального радиуса в угол графика
text(1.1 * max(abs(x_opt(1:2:end))), 1.1 * max(abs(x_opt(2:2:end))), ['R = ', num2str(x_opt(end))], 'Color', 'black', 'FontSize', 12, 'FontWeight', 'bold');

axis equal;
hold off;

pause;
