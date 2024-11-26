u = imread('eye.png');
u = double(u) / 255;
figure(1); clf;
imagesc(u);
%pcolor(u);
caxis([0 1]);

colormap(gray);
axis equal;
title('Initial picture');
axis tight;
exportgraphics(gcf, "eye_unblurred.png", 'Resolution', 300);

[n, n2] = size(u);
if (n ~= n2)
    error('by default, this only supports square images');
end

u = reshape(u, n*n, 1); 
u = u0

h = 1;  
k = 0.1; 
s = (h:h:n)';  

[xx, yy] = meshgrid(s, s);
x = xx(:); 
y = yy(:);

N = length(s);
e = ones(N, 1);
L1d = spdiags([e -2 * e e], [-1 0 1], N, N);
L1d(1, 1) = -1; 
L1d(1, 2) = 1;
L1d(N, N) = -1;  
L1d(N, N-1) = 1;
L1d = 1 / h^2 * L1d; 
I1d = speye(size(L1d));

L = kron(I1d, L1d) + kron(L1d, I1d);  

numsteps = 10;  
t = 0;


for n = 1:numsteps
    uu = reshape(u, size(xx));  
    
    u = u + k * (L * u); 
    t = t + k;  
end

u2 = reshape(u, n2, n2); 
figure(2); clf;
imagesc(u2);
%pcolor(u2);
caxis([0 1]);
colormap(gray);
axis equal;
axis tight;
title('Final picture');
exportgraphics(gcf, "eye_blurred.png", 'Resolution', 300);