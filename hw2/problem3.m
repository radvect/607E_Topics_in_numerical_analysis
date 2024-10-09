x = linspace(0, 8, 1000); 

L01 = (x-1)./(0-1)*2 + (x-0)./(1-0)*7;
L12 = (x-3)./(1-3)*7 + (x-1)./(3-1)*5;
L23 = (x-6)./(3-6)*5 + (x-3)./(6-3)*0;

L02 = (x - 3) ./ (0 - 3) .* L01 + (x - 0) ./ (3 - 0) .* L12;
L13 = (x - 6) ./ (1 - 6) .* L12 + (x - 1) ./ (6 - 1) .* L23;

L03 = (x - 6) ./ (0 - 6) .* L02 + (x - 0) ./ (6 - 0) .* L13;

figure;
hold on;


plot(x, L01, 'b-', 'LineWidth', 1.5, 'DisplayName', 'L01');
plot(x, L12, 'g--', 'LineWidth', 1.5, 'DisplayName', 'L12');
plot(x, L23, 'r-.', 'LineWidth', 1.5, 'DisplayName', 'L23');


plot(x, L02, 'm:', 'LineWidth', 2, 'DisplayName', 'L02');
plot(x, L13, 'c-', 'LineWidth', 2, 'DisplayName', 'L13');
plot(x, L03, 'k-', 'LineWidth', 2, 'DisplayName', 'L03');


plot(0, 2, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'HandleVisibility', 'off');
plot(1, 7, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'HandleVisibility', 'off');
plot(3, 5, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'HandleVisibility', 'off');
plot(6, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'HandleVisibility', 'off');


xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y', 'FontSize', 12, 'FontWeight', 'bold');
legend('show', 'Location', 'best', 'FontSize', 10);
grid on;

hold off;
saveas(gcf, 'C:\\Users\\Pavel\\Documents\\MATLAB\\LinearPolynomsComparison.png'); 




figure;
hold on;

plot(x, L03, 'k-', 'LineWidth', 2, 'DisplayName', 'L03');

plot(0, 2, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'HandleVisibility', 'off');
plot(1, 7, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'HandleVisibility', 'off');
plot(3, 5, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'HandleVisibility', 'off');
plot(6, 0, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'HandleVisibility', 'off');

xlabel('x', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('y', 'FontSize', 12, 'FontWeight', 'bold');
legend('show', 'Location', 'best', 'FontSize', 10);
grid on;

hold off;
saveas(gcf, 'C:\\Users\\Pavel\\Documents\\MATLAB\\final_L03.png'); 

