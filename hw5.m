clear all; close all; clc;
%%
X = VideoReader('badminton.mp4');
video = [];
while hasFrame(X)
    data = readFrame(X);
    data = rgb2gray(data);
    data = data(:);
    video = [video data];
end
video = double(video);

width = X.width;
height = X.height;
dt = 1 / X.FrameRate;

X1 = video(:, 1:end-1); X2 = video(:, 2:end);
frame = size(video, 2) - 1;
[U, S, V] = svd(X1, 'econ');
%plot(diag(S)/sum(diag(S)), 'ro');
%%
r = 1;
Ur = U(:, 1:r);
Sr = S(1:r, 1:r);
Vr = V(:, 1:r);

Atilde = Ur'*X2*Vr / Sr;
[W, D] = eig(Atilde);
Phi = X2*Vr/ Sr*W;
lambda = diag(D);
omega = log(lambda) / dt;
[val, ind] = sort(abs(omega));

omega = omega(ind(1));
%%
t = linspace(0, X.Duration, size(video, 2));
x1 = video(:, 1);
b = Phi\x1;
time_dynamics = zeros(r,length(t));
for iter = 1:length(t)
    temp = (b.*(exp(omega*t(iter))).');
    time_dynamics(:,iter) = sum(temp, 2);
end

X_lowRank = Phi*time_dynamics;
X_sparse = video - abs(X_lowRank);

R = zeros(length(X_sparse), size(X_sparse, 2));

% 
% for i = 1:length(X_sparse)
%     for j = 1:frame
%         if X_sparse(i, j) < 0
%             R(i,j) = X_sparse(i,j);
%         end
%     end
% end
X_sparse = X_sparse - R;
X_lowRank = abs(X_lowRank) + R;

video_sparse = reshape(X_sparse, [height, width, frame+1]);
video_lowRank = reshape(X_lowRank, [height, width, frame+1]);
%% make video
videoPlayer = vision.VideoPlayer;
for j = 1:frame
    videoFrame = video_sparse(:,:,j);
    %videoFrame = video_lowRank(:,:,j);
    videoPlayer(uint8(videoFrame));
    pause(dt);
end

%% plot results
i = 100;
figure(1)
subplot(3,3,1);
pcolor(flip(reshape(video(:,i), [height, width]))); shading interp; colormap(gray);axis off;
title('Original film at frame = 100');
subplot(3,3,2);
pcolor(flip(video_lowRank(:,:,i))); shading interp; colormap(gray);axis off;
title('Background at frame = 100');
subplot(3,3,3);
pcolor(flip(video_sparse(:,:,i))); shading interp; colormap(gray);axis off;
title('Foreground at frame = 100');
%%
%i = 150;
subplot(3,3,4);
pcolor(flip(reshape(video(:,i), [height, width]))); shading interp; colormap(gray);axis off;
title('Original film at frame = 100');
subplot(3,3,5);
pcolor(flip(video_lowRank(:,:,i))); shading interp; colormap(gray);axis off;
title('Background at frame = 100');
subplot(3,3,6);
pcolor(flip(video_sparse(:,:,i))); shading interp; colormap(gray);axis off;
title('Foreground at frame = 100');
%%
%i = 238;
subplot(3,3,7);
pcolor(flip(reshape(video(:,i), [height, width]))); shading interp; colormap(gray);axis off;
title('Original film at frame = 100');
subplot(3,3,8);
pcolor(flip(video_lowRank(:,:,i))); shading interp; colormap(gray);axis off;
title('Background at frame = 100');
subplot(3,3,9);
pcolor(flip(video_sparse(:,:,i))); shading interp; colormap(gray);axis off;
title('Foreground at frame = 100');


%%
i = 200;
subplot(1,2,1);
pcolor(flip(video_sparse(:,:,i))); shading interp; colormap(gray);axis off;
title({['Foreground at frame = 200'];['before subtracting R']});


for i = 1:length(X_sparse)
    for j = 1:330
        if X_sparse(i, j) < 0
            R(i,j) = X_sparse(i,j);
        end
    end
end
X_sparse = X_sparse - R;
X_lowRank = abs(X_lowRank) + R;

video_sparse = reshape(X_sparse, [height, width, frame+1]);
video_lowRank = reshape(X_lowRank, [height, width, frame+1]);

i = 200;
subplot(1,2,2);
pcolor(flip(video_sparse(:,:,i))); shading interp; colormap(gray);axis off;
title({['Foreground at frame = 200'];['after subtracting R']});
