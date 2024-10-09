n=0;

function optimize_circles()
    % Количество малых кругов
    n = 6; % Замените на нужное количество кругов

    % Радиусы малых кругов (все равны 1)
    radii = ones(1, n);

    % Начальные приближения: размещаем круги по окружности внутри контейнера
    theta = linspace(0, 2*pi, n+1);
    theta(end) = [];
    initial_R = 2;  % Начальный радиус контейнера (должен быть >= 1)
    initial_vars = zeros(1, 2*n + 1);
    for i = 1:n
        initial_vars(2*i - 1) = (initial_R - 1) * cos(theta(i));
        initial_vars(2*i) = (initial_R - 1) * sin(theta(i));
    end
    initial_vars(end) = initial_R;

    % Запуск оптимизации
    options = optimset('MaxIter', 10000, 'MaxFunEvals', 100000, 'TolFun', 1e-8, 'TolX', 1e-8,'Display', 'iter');
    [opt_vars, opt_val] = fminsearch(@(vars) objective_function(vars, radii), initial_vars, options);

    % Результаты
    x_coords = opt_vars(1:2:2*n - 1);
    y_coords = opt_vars(2:2:2*n);
    R_opt = opt_vars(end);

    % Вывод результатов
    fprintf('Оптимальный радиус контейнера: %f\n', R_opt);
    for i = 1:n
        fprintf('Круг %d: x = %f, y = %f, r = %f\n', i, x_coords(i), y_coords(i), radii(i));
    end

    % Визуализация
    plot_circles(x_coords, y_coords, radii, R_opt);
end

function total_cost = objective_function(vars, radii)
    n = length(radii);
    P = 1e6;  % Коэффициент штрафа

    x = vars(1:2:2*n - 1);
    y = vars(2:2:2*n);
    R = vars(end);

    % Целевая функция: минимизировать R
    total_cost = R;

    % Штраф за пересечение кругов
    penalty_overlap = 0;
    for i = 1:n
        for j = i+1:n
            dist_sq = (x(i) - x(j))^2 + (y(i) - y(j))^2;
            min_dist_sq = (radii(i) + radii(j))^2;
            overlap = max(0, min_dist_sq - dist_sq);
            penalty_overlap = penalty_overlap + overlap * P;
        end
    end

    % Штраф за выход за пределы контейнера
    penalty_outside = 0;
    for i = 1:n
        dist_to_center = sqrt(x(i)^2 + y(i)^2) + radii(i);
        outside = max(0, dist_to_center - R);
        penalty_outside = penalty_outside + outside * P;
    end

    % Суммарный штраф
    total_cost = total_cost + penalty_overlap + penalty_outside;
end

function plot_circles(x, y, radii, R)
    figure;
    hold on;
    axis equal;

    % Рисуем контейнерный круг
    theta_plot = linspace(0, 2*pi, 100);
    plot(R * cos(theta_plot), R * sin(theta_plot), 'b--', 'LineWidth', 2);

    % Рисуем малые круги
    for i = 1:length(radii)
        circle_x = x(i) + radii(i) * cos(theta_plot);
        circle_y = y(i) + radii(i) * sin(theta_plot);
        fill(circle_x, circle_y, 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    end

    xlabel('X');
    ylabel('Y');
    title('Упаковка кругов в круге');
    grid on;
    hold off;
    pause
end
optimize_circles();

