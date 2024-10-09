pkg load optim; % Загрузка пакета для оптимизации

function d = circle_distance(x1, y1, x2, y2)
  % Функция для расчета расстояния между центрами двух кругов
  d = sqrt((x2 - x1)^2 + (y2 - y1)^2);
end

function penalty = packing_cost(positions, r_small, R_big, N)
  % Функция стоимости для оптимизации (наказание за пересечение или выход за пределы)
  penalty = 0;

  % Проверка, что каждый круг внутри большого круга
  for i = 1:N
    x = positions(2*i - 1);
    y = positions(2*i);
    if circle_distance(0, 0, x, y) + r_small > R_big
      penalty += 1000; % Наказание за выход за пределы
    end
  end

  % Проверка на пересечения между кругами
  for i = 1:N
    for j = i+1:N
      x1 = positions(2*i - 1);
      y1 = positions(2*i);
      x2 = positions(2*j - 1);
      y2 = positions(2*j);
      if circle_distance(x1, y1, x2, y2) < 2 * r_small
        penalty += 1000; % Наказание за пересечение
      end
    end
  end
end

function packed_positions = pack_circles(N, r_small, R_big)
  % Функция для упаковки N кругов радиуса r_small в круг радиуса R_big
  % Инициализация случайных начальных позиций
  initial_positions = rand(1, 2 * N) * 2 * R_big - R_big;

  % Ограничения на границы (в пределах большого круга)
  lb = -R_big * ones(1, 2 * N);
  ub = R_big * ones(1, 2 * N);

  % Использование оптимизатора для минимизации штрафа
  options = optimset('Display', 'iter');
  packed_positions = fmincon(@(positions) packing_cost(positions, r_small, R_big, N), initial_positions, [], [], [], [], lb, ub, [], options);
end

% Пример использования:
N = 10;        % Количество маленьких кругов
r_small = 1;   % Радиус маленького круга
R_big = 10;    % Радиус большого круга

% Запуск упаковки
positions = pack_circles(N, r_small, R_big);

% Отрисовка результата
theta = linspace(0, 2*pi, 100);

hold on;
% Рисуем большой круг
plot(R_big * cos(theta), R_big * sin(theta), 'b');

% Рисуем маленькие круги
for i = 1:N
  x = positions(2*i - 1);
  y = positions(2*i);
  plot(x + r_small * cos(theta), y + r_small * sin(theta), 'r');
end

axis equal;
hold off;
