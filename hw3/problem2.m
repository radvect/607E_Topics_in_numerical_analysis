format long 

function w = calculate_w(xk)
    n = length(xk); 
    w = zeros(1, n);

    for k = 1:n
        B_xk = 1; 
        for i = 1:n
            if i ~= k
                B_xk = B_xk * (xk(k) - xk(i)); 
            end
        end
        
        w(k) = 1 / B_xk; 
    end
end



function beta_vector = bary(xk, x)
    w = calculate_w(xk)
    idx = find(abs(x - xk) < eps);
    n = length(w)
    beta_vector = zeros(1, n)
    if ~isempty(idx)
        beta_vector(idx) = 1;
    else 
        denominator = 0
        for i = 1:n
            denominator = denominator + w(i)/(x - xk(i))
        end
        for i = 1:n
            beta_vector(i) = w(i)/(x-xk(i))/denominator
        end
    end
end


xk = [0 1 3]
x = 0
beta = bary(xk, x)
f = [15 8 0]
finterp = sum(beta .* f)


