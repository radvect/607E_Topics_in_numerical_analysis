format long

h_array = [];
error_array = [];
fpp_theory = sin(1.2)
for k = 0:.5:8
    h = 10^(-k)
    x = 1.2;
    fpp_exp = 1/h/h/h/h*(sin(x-2*h)-4*sin(x-h)+6*sin(x)-4*sin(x+h)+sin(x+2*h))
    error_array = [error_array, abs(fpp_theory - fpp_exp)]
    h_array = [h_array, h]
end
loglog(h_array, error_array)
hold on
title(sprintf("error(h)"))
xlabel('h')
ylabel('error')
hold off
saveas(gcf, 'C:\\Users\\Pavel\\Documents\\MATLAB\\hw3\\pic\\error_fourth.png');