x = 0

function isInside = isPointInside(u, v, r, x, y)
    distance = sqrt((x - u)^2 + (v - y)^2);
    if distance < r
        isInside = -1;
    elseif abs(distance-r)<1e-10
        isInside = 0;
    else
        isInside = 1;
    end
end


function distance = distance_two_points(x1, x2, y1, y2)
    distance = sqrt((x1 - x2)^2 + (y1 - y2)^2);
end

function distance_value= sdcircle_1(u, v, r, x, y)

    dphi = 0.1;
    distance_value = inf;
    distance_value_new = 0
    while abs((distance_value) - (distance_value_new)) > 1e-3
        distance_value = distance_value_new;
        phi = 0:dphi:2*pi;
        x_circle = u + r * cos(phi);

        y_circle = v + r * sin(phi);
        distance_new = zeros(length(x_circle), 1);
        for i = 1:length(x_circle)
            distance_new(i) = distance_two_points(x, x_circle(i), y, y_circle(i));
        end
        [distance_value_new, min_index] = min(distance_new)
        dphi = dphi/10.;
    end
    distance_value = distance_value* isPointInside(u, v, r, x, y);

end

function distance_value= sdcircle_2(u, v, r, x, y)

    distance_value = abs((sqrt((u - x)^2 + (v - y)^2))-r)* isPointInside(u, v, r, x, y);
end



function [X, Y]  = compare_contour_maps()
    % Параметры окружностей
    u = [2, -2, 0];
    v = [2, 2, -1];
    r = [1.25, 1.5, 3];

    % Создаем сетку точек
    [X, Y] = meshgrid(linspace(-5, 5, 50), linspace(-5, 5, 50));

    % Первая часть: построение для sdcircle_1
    subplot(1, 2, 1);
    f_min_1 = inf * ones(size(X));
    for i = 1:3
        f_circle_1 = arrayfun(@(x, y) sdcircle_1(u(i), v(i), r(i), x, y), X, Y);
        f_min_1 = min(f_min_1, f_circle_1);
    end
    contourf(X, Y, f_min_1, 30); % Контуры с заполнением
    hold on;

    % Построение окружностей
    theta = linspace(0, 2*pi, 100);
    for i = 1:3
        x_circle = u(i) + r(i) * cos(theta);
        y_circle = v(i) + r(i) * sin(theta);
        plot(x_circle, y_circle, 'k--', 'LineWidth', 2); % Строим окружности как пунктирные линии
    end
    colorbar;
    xlabel('x');
    ylabel('y');
    title("Distance for sdcircle_1");
    axis equal;

    % Вторая часть: построение для sdcircle_2
    subplot(1, 2, 2);
    f_min_2 = inf * ones(size(X));
    for i = 1:3
        f_circle_2 = arrayfun(@(x, y) sdcircle_2(u(i), v(i), r(i), x, y), X, Y);
        f_min_2 = min(f_min_2, f_circle_2);
    end
    contourf(X, Y, f_min_2, 30); % Контуры с заполнением
    hold on;

    % Построение окружностей для второго графика
    for i = 1:3
        x_circle = u(i) + r(i) * cos(theta);
        y_circle = v(i) + r(i) * sin(theta);
        plot(x_circle, y_circle, 'k--', 'LineWidth', 2); % Строим окружности как пунктирные линии
    end
    colorbar;
    xlabel('x');
    ylabel('y');
    title("Distance for sdcircle_2");
    axis equal;

    pause;
end

[X,Y] = compare_contour_maps();


%[X,Y] = compare_contour_maps()

