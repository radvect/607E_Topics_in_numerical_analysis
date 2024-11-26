u = imread('problems.png');
if ndims(u) == 3
    u = rgb2gray(u); 
end

u = double(u) / 255;


figure(1); clf;
imagesc(u);
%pcolor(u);
caxis([0 1]);

colormap(gray);
axis equal;
title('Initial picture');
axis tight;
exportgraphics(gcf, "IMAGE_book.png", 'Resolution', 300);

[n, n2] = size(u);
if (n ~= n2)
    error('by default, this only supports square images');
end

u0 = reshape(u, n*n, 1); 
disp(max(u0)*255);
disp(min(u0)*255);
u = u0;
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

numsteps = 100;  
t = 0;


for n = 1:numsteps
    uu = reshape(u, size(xx));  
    
    u = u + k * (L * u); 
    t = t + k;  
end

edge_map = u0 - u;

acute = u0 + edge_map;
disp(max(acute)*255);
disp(min(acute)*255);
u2 = reshape(acute, n2, n2); 

figure(2); clf;
imagesc(u2);
%pcolor(u2);
caxis([0 1]);
colormap(gray);

axis equal;
axis tight;
title('After applying mask');
exportgraphics(gcf, "book10.png", 'Resolution', 300);