
function y = next_single(x)
    y = single(x + eps(x));
end


x = single(16777208);

for i = 1:20
    binary_repr = dec2bin(typecast(x, 'uint32'), 32); 
    fprintf('%s\t\t%.7f\n', binary_repr, x);
    x = next_single(x);
end