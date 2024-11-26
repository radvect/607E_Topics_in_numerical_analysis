format long
N = 65536;
SIN_TABLE = zeros(1, N, 'single');
for i = 1:N
    sin_value = sin((i-1) * (2*pi) / N);
    SIN_TABLE(i) = single(sin_value); 
end

function y = next_single(x)
    y = single(x + eps(x));
end

yaw = single(5850000);
max_var = -1;
best_yaw = yaw;

for i = 1:1000000
    const1 = single(0.017453292);
    const2 = single(10430.378);
    value = yaw * const1;

    sindex = rem(int64(value * const2), N) + 1;
    cindex = rem(int64(value * const2 + 16384), N) + 1;
    s = SIN_TABLE(sindex);
    c = SIN_TABLE(cindex);

    new_var = c*c + s*s;
    
    if new_var > max_var
        fprintf('value*const2: %.9sf\n', (value * const2));
        fprintf('value*const2+16384: %.9sf\n', (value * const2)+16834);

        int64(value * const2)+16384

        int64(value * const2+16384)
        max_var = new_var;
        best_yaw = yaw;
    end

    yaw = next_single(yaw);
end

fprintf('Best yaw: %.10f\n', best_yaw);
fprintf('Maximum c^2 + s^2: %.10f\n', max_var);

