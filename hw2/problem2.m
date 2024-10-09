x = linspace(-10, 10, 1000); 

y0 = ((x-1).*(x-3).*(x-6)) / -18;
y1 = (x.*(x-3).*(x-6)) / 10;
y2 = (x.*(x-1).*(x-6)) / -18;
y3 = ((x).*(x-1).*(x-3)) / 90;

figure;
hold on;

plot(x, y0, 'r-', 'DisplayName', 'L(4,0)');
plot(x, y1, 'g--', 'DisplayName', 'L(4,1)');
plot(x, y2, 'b-.', 'DisplayName', 'L(4,2)');
plot(x, y3, 'm:', 'DisplayName', 'L(4,3)'); 


xlabel('x');
ylabel('y');
legend('show');
grid on;

hold off;
saveas(gcf, 'C:\\Users\\Pavel\\Documents\\MATLAB\\polynomials_plot.png'); 