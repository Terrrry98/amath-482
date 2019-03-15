clear all; close all; clc;
load Testdata
L=15; % spatial domain
n=64; % Fourier modes
x2=linspace(-L,L,n+1); x=x2(1:n); y=x; z=x;
k=(2*pi/(2*L))*[0:(n/2-1) -n/2:-1]; ks=fftshift(k);
[X,Y,Z]=meshgrid(x,y,z);
[Kx,Ky,Kz]=meshgrid(ks,ks,ks);

%% 1.

uave = zeros(n,n,n);
for i = 1:20
    Un(:,:,:)=reshape(Undata(i,:),n,n,n);
    Uf = fftn(Un);
    uave = uave + Uf;
end

% Normalize
uave = abs(fftshift(uave))./ max(max(max(abs(uave))));

% Find the index of the center frequency
[C, I] = max(uave(:));
[X_ind, Y_ind, Z_ind] = ind2sub([n,n,n], I);
% X_ind = 28; Y_ind = 42; Z_ind = 33;

% Find the corresponding frequency
x_freq = Kx(X_ind, Y_ind, Z_ind);
y_freq = Ky(X_ind, Y_ind, Z_ind);
z_freq = Kz(X_ind, Y_ind, Z_ind);

% Plot isosurface
figure(1)
isosurface(X, Y, Z, abs(uave), 0.4); grid on;
axis([-20 20 -20 20 -20 20]);


%% 2.

% Create a filter
filter = exp(-0.2.*(Kx - x_freq).^2 - 0.2.*(Ky - y_freq).^2 - 0.2.*(Kz - z_freq).^2);
filter = fftshift(filter);

% Create an empty matrix to store the path
path = zeros(20,3);

for i = 1:20
    Un(:,:,:)=reshape(Undata(i,:),n,n,n);
    utn = fftn(Un);
    utnf = utn .* filter;
    unf = ifftn(utnf);
    
    % Find x, y, z coordinates
    [M, I_2] = max(unf(:));
    [X_path, Y_path, Z_path] = ind2sub([n,n,n], I_2);
    x_path = X(X_path, Y_path, Z_path);
    y_path = Y(X_path, Y_path, Z_path);
    z_path = Z(X_path, Y_path, Z_path);
    
    path(i, 1) = x_path;
    path(i, 2) = y_path;
    path(i, 3) = z_path;
end

% Plot the path
plot3(path(:,1), path(:,2), path(:,3), 'LineWidth', 2); grid on;
title('Path of the marble');
xlabel('x'); ylabel('y'); zlabel('z');
hold on;
%% 3.
% 20th data
plot3(path(20,1), path(20,2), path(20,3), '*', 'MarkerSize', 10);

hold off;