format long

u = 2%+eps*10e5;
v = 1%+eps*10e5;
r = 4%+eps*10e5;
theta = pi/6 %0%pi/3% pi/6+eps*10e5;

x = -5:1:5
y = -5:1:5


function [P1_rot,P2_rot,P3_rot,P4_rot] = get_rotated_vertex_coordinates(u, v, r, theta)
    P1 = [u - r/2, v - r/2];
    P2 = [u + r/2, v - r/2];
    P3 = [u + r/2, v + r/2];
    P4 = [u - r/2, v + r/2];
    R = [cos(theta), -sin(theta); sin(theta), cos(theta)];

    rotate_point = @(P) R * (P' - [u; v]) + [u; v];


    P1_rot = rotate_point(P1);
    P2_rot = rotate_point(P2);
    P3_rot = rotate_point(P3);
    P4_rot = rotate_point(P4);

end


function isInside = isPointInside(u, v, r, theta, x, y)
    point = [x; y];
    R_inv = [cos(-theta), -sin(-theta); sin(-theta), cos(-theta)];
    point_rel = point - [u; v];
    point_rot = R_inv * point_rel;

    if abs(point_rot(1)) < r/2 && abs(point_rot(2)) < r/2
        isInside = -1;
    elseif abs(point_rot(1)) == r/2 || abs(point_rot(2)) == r/2
        isInside = 0;
    else
        isInside=1
    end
end
function distance = distance_two_points(x1, y1, x2, y2)
    distance = sqrt((x1 - x2)^2 + (y1 - y2)^2);
end

function distance_value = sdsquare_add(u, v, r, theta, x, y)
    npoints = 2;
    distance_value = inf;
    distance_value_new = 0;

    [P1_rot, P2_rot, P3_rot, P4_rot] =  get_rotated_vertex_coordinates(u, v, r, theta);
    %figure;
    %hold on;
    %plot([P1_rot(1), P2_rot(1), P3_rot(1), P4_rot(1), P1_rot(1)], [P1_rot(2), P2_rot(2), P3_rot(2), P4_rot(2), P1_rot(2)], 'r-', 'LineWidth', 2);
    %plot(x, y, 'bo-')
    %pause;
    while abs((distance_value) - (distance_value_new)) > 1e-3
        distance_value = distance_value_new

        x_bottom = linspace(P1_rot(1), P2_rot(1), npoints);
        y_bottom = linspace(P1_rot(2), P2_rot(2), npoints);

        x_right = linspace(P2_rot(1), P3_rot(1), npoints);
        y_right = linspace(P2_rot(2), P3_rot(2), npoints);

        x_top = linspace(P3_rot(1), P4_rot(1), npoints);
        y_top = linspace(P3_rot(2), P4_rot(2), npoints);

        x_left = linspace(P4_rot(1), P1_rot(1), npoints);
        y_left = linspace(P4_rot(2), P1_rot(2), npoints);

        points_x = [x_bottom, x_right, x_top,  x_left ];
        points_y =[y_bottom, y_right, y_top,  y_left ];

        distance_new = zeros(length(points_x), 1);
        for i = 1:length(points_x)
            distance_new(i) = distance_two_points(x, y, points_x(i),  points_y(i));
        end

        [distance_value_new, min_index] = min(distance_new);

        npoints = npoints*2
        abs((distance_value) - (distance_value_new));
    end
    distance_value = distance_value* isPointInside(u, v, r, theta, x, y);

end


function distance_value = sdsquare(u, v, r, theta, x, y)
    point = [x,y];
    R_inv = [cos(-theta), -sin(-theta); sin(-theta), cos(-theta)];

    point_rot = R_inv * (point' - [u; v]) + [u; v];
    [P1_rot,P2_rot,P3_rot,P4_rot] = get_rotated_vertex_coordinates(u, v, r,0)


    if (point_rot(1) > min([P1_rot(1), P2_rot(1), P3_rot(1), P4_rot(1)]) && point_rot(1) < max([P1_rot(1), P2_rot(1), P3_rot(1), P4_rot(1)]))
        in_x = true
    else
        in_x = false
    end
    if (point_rot(2) > min([P1_rot(2), P2_rot(2), P3_rot(2), P4_rot(2)]) && point_rot(2) < max([P1_rot(2), P2_rot(2), P3_rot(2), P4_rot(2)]))
        in_y = true
    else
        in_y = false
    end


    if(~ in_x && ~ in_y)
        distance_value = min([distance_two_points(point_rot(1),point_rot(2),P1_rot(1), P1_rot(2)), distance_two_points(point_rot(1),point_rot(2),P2_rot(1), P2_rot(2)),distance_two_points(point_rot(1),point_rot(2),P3_rot(1), P3_rot(2)),distance_two_points(point_rot(1),point_rot(2),P4_rot(1), P4_rot(2))])
    elseif(~ in_x && in_y)
        distance_value = min([distance_two_points(point_rot(1),point_rot(2),P1_rot(1), point_rot(2)), distance_two_points(point_rot(1),point_rot(2),P2_rot(1), point_rot(2)),distance_two_points(point_rot(1),point_rot(2),P3_rot(1), point_rot(2)),distance_two_points(point_rot(1),point_rot(2),P4_rot(1), point_rot(2))])
    elseif(in_x && ~ in_y)
     distance_value = min([distance_two_points(point_rot(1),point_rot(2),point_rot(1), P1_rot(2)), distance_two_points(point_rot(1),point_rot(2),point_rot(1), P2_rot(2)),distance_two_points(point_rot(1),point_rot(2),point_rot(1), P3_rot(2)),distance_two_points(point_rot(1),point_rot(2),point_rot(1), P4_rot(2))])
    elseif(in_x && in_y)
      x_min_distance = min([distance_two_points(point_rot(1),point_rot(2),P1_rot(1), point_rot(2)), distance_two_points(point_rot(1),point_rot(2),P2_rot(1), point_rot(2)),distance_two_points(point_rot(1),point_rot(2),P3_rot(1), point_rot(2)),distance_two_points(point_rot(1),point_rot(2),P4_rot(1), point_rot(2))])
      y_min_distance = min([distance_two_points(point_rot(1),point_rot(2),point_rot(1), P1_rot(2)), distance_two_points(point_rot(1),point_rot(2),point_rot(1), P2_rot(2)),distance_two_points(point_rot(1),point_rot(2),point_rot(1), P3_rot(2)),distance_two_points(point_rot(1),point_rot(2),point_rot(1), P4_rot(2))])
      distance_value = -min([x_min_distance, y_min_distance])
  else
    distance_value = 0
    end

  if (isPointInside(u, v, r, 0, point_rot(1), point_rot(2))==0)
      distance_value=0
  endif


end

[x, y] = meshgrid(-5:0.1:5, -5:0.1:5);


Z = zeros(size(x));


for i = 1:size(x, 1)
     for j = 1:size(x, 2)
        Z(i, j) = sdsquare(u, v, r, theta , x(i, j), y(i, j));
     end
end


figure;
contourf(x, y, Z, 20);
colorbar;
title('Distance to the square boundary');
xlabel('x');
ylabel('y');
pause;
