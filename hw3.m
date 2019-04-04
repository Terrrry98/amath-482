%% normal
load cam1_1;
load cam2_1;
load cam3_1;
load cam1_2;
load cam2_2;
load cam3_2;
load cam1_3;
load cam2_3;
load cam3_3;
load cam1_4;
load cam2_4;
load cam3_4;
%%
numFrames1 = size(vidFrames1_1, 4);
numFrames2 = size(vidFrames2_1, 4);
numFrames3 = size(vidFrames3_1, 4);
%%
y_1 = [];
x_1 = [];
for j=1:numFrames1
    X1 = double(rgb2gray(vidFrames1_1(:,:,:,j)));
    X1(:, 1:300) = 0; X1(:, 400:end) = 0; X1(1:200, :) = 0; X1(400:end, :)=0;
    [V, I] = max(X1(:));
    [aaa, bbb] = find(X1 >= V * 11 / 12);
    y_1(j) = mean(aaa);
    x_1(j) = mean(bbb);
end
[v, i] = max(y_1(1:50));
x_1 = x_1(30:end);
y_1 = y_1(30:end);
%%
y_2 = [];
x_2 = [];
for j=1:numFrames2
    X2 = double(rgb2gray(vidFrames2_1(:,:,:,j)));
    X2(:, 1:220) = 0; X2(:, 350:end) = 0; X2(1:100, :) = 0; X2(350:end, :)=0;
    [V, I] = max(X2(:));
    [aaa, bbb] = find(X2 >= V * 11 / 12);
    y_2(j) = mean(aaa);
    x_2(j) = mean(bbb);
end
[v, i] = max(y_2(1:41));
x_2 = x_2(i:end);
y_2 = y_2(i:end);
%%
x_3 = [];
y_3 = [];
for j=1:numFrames3
    X3 = double(rgb2gray(vidFrames3_1(:,:,:,j)));
    X3(1:230, :) = 0; X3(350:end, :) = 0; X3(:, 1:200) = 0; X3(:, 480:end)=0;
    [V, I] = max(X3(:));
    [aaa, bbb] = find(X3 >= 11 / 12 * V);
    y_3(j) = mean(aaa);
    x_3(j) = mean(bbb);
end

[v, i] = max(y_3(1:41));
x_3 = x_3(i:end);
y_3 = y_3(i:end);
%%
l = min([length(y_1), length(y_2), length(x_3)]);
X = [x_1(1:l); y_1(1:l); x_2(1:l); y_2(1:l); x_3(1:l); y_3(1:l)];
[m, n] = size(X);
mn = mean(X, 2);
X=X-repmat(mn,1,n); % subtract mean --- put in the same scale
[u, s, v] = svd(X, 'econ');
figure(1)
plot(diag(s)./sum(diag(s)), 'ro')
xlabel('Principal component'); ylabel('energy percentage(%)');
v = v*s;
figure(2)
plot(v(:,1)); hold on;
plot(v(:,2));
plot(v(:,3));
plot(v(:,4));
plot(v(:,5)); 
plot(v(:,6)); hold off;
xlabel('time(frame numbers)')
ylabel('position')
title('Case 1')
legend('component1', 'component2', 'component3', 'component4', 'component5', 'component6')