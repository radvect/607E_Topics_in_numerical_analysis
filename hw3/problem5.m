k = 17.0652165601579625588917206249 / 29000;  
mu = 0.012277471;      
mu_bar = 1 - mu;

u11 = 0.994; 
u12 = 0.0;
u21 = 0.0; 
u22 = -2.00158510637908252240537862224;
U_n = [u11; u12; u21; u22];

num_iter = 29000;

u11_values = zeros(num_iter, 1);
u12_values = zeros(num_iter, 1);

for i = 1:num_iter

    u11_values(i) = U_n(1);
    u12_values(i) = U_n(2);

 
    D1 = ((U_n(1) + mu)^2 + (U_n(2))^2)^(3/2);
    D2 = ((U_n(1) - mu_bar)^2 + (U_n(2))^2)^(3/2);


    A = [1,  0,  k,  0;
         0,  1,  0,  k;
         k*(1 - mu_bar/D1 - mu/D2),  0,  1,  2*k;
         0,  k*(1 - mu_bar/D1 - mu/D2),  -2*k,  1];

    B = [0; 0; k*mu*mu_bar*(-1/D1 + 1/D2); 0];


    U_next = A * U_n + B;

    U_n = U_next;

end


figure;
plot(u11_values, u12_values,'b', 'LineWidth', 1.5);
xlabel('y1');
ylabel('y2');
title('Trajectory of the spacecraft');
grid on;
saveas(gcf, 'C:\\Users\\Pavel\\Documents\\MATLAB\\hw3\\pic\\spacecraft1.png');