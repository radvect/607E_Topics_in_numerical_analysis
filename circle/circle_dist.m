u = 3/2+eps*10e5;
v = pi+eps*10e5;
r = e+eps*10e5;
x = 10.+eps*10e5;
y = 1.+eps*10e5;
eps_value = eps;


function isInside = isPointInside(u, v, r, x, y)
    distance = sqrt((x - u)^2 + (v - y)^2);
    if distance < r
        isInside = 1;
    elseif abs(distance-r)<1e-10
        isInside = 0;
    else
        isInside = -1;
    end
end


function distance = distance_two_points(x1, x2, y1, y2)
    distance = sqrt((x1 - x2)^2 + (y1 - y2)^2);
end

function distance_value= sdcircle(u, v, r, x, y)

    dphi = 1.;
    distance_value = inf;
    distance_value_new = 0
    while abs((distance_value) - (distance_value_new)) > 1e-10
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

result = sdcircle(u, v, r, x, y);
disp(result);
