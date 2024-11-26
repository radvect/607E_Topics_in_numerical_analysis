format long
N = 65536;
SIN_TABLE = zeros(1, N, 'single');
for i=1:1:N
    sin_value = sin((i-1) * (2*pi) / N);
    SIN_TABLE(i) = single(sin_value); 
end
fprintf('value_in_computations_type:%s\n',  class(sin_value)); 
fprintf('SIN_TABLE(1):%.10f,type:%s\n', SIN_TABLE(1), class(SIN_TABLE(1))); 
fprintf('SIN_TABLE(N):%.10f,type:%s\n', SIN_TABLE(N), class(SIN_TABLE(N))); 
fprintf('SIN_TABLE(43):%.10f,type:%s\n', SIN_TABLE(43), class(SIN_TABLE(43))); 