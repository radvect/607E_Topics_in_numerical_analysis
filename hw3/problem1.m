x = linspace(-3, 3, 500);  
y = linspace(-3, 3, 500);  
[X, Y] = meshgrid(x, y);  


klambda = X + 1i * Y;


F1 = abs(klambda + sqrt(klambda.^2 + 1));  
F2 = abs(klambda - sqrt(klambda.^2 + 1));  


figure;


subplot(1, 2, 1);
contourf(X, Y, F1 <= 1, [1 1], 'LineWidth', 1.5);
colormap('winter');  
axis equal;
title('|k\lambda + \surd(k\lambda^2 + 1)| \leq 1');
xlabel('Re(k\lambda)');
ylabel('Im(k\lambda)');
grid on;


subplot(1, 2, 2);
contourf(X, Y, F2 <= 1, [1 1], 'LineWidth', 1.5);
colormap('autumn');  
axis equal;
title('|k\lambda - \surd(k\lambda^2 + 1)| \leq 1');
xlabel('Re(k\lambda)');
ylabel('Im(k\lambda)');
grid on;



saveas(gcf, 'C:\\Users\\Pavel\\Documents\\MATLAB\\hw3\\pic\\reim.png');
