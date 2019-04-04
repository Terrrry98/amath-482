%% load data.
path = './song/test2';
folder = dir(path);

A1 = [];
for i = 3:length(folder)
    filename = folder(i).name;
    [song, Fs1] = audioread(strcat(path, '/', filename)); %load data
    song = song'/ 2;
    song = song(1,:) + song(2,:);   % merge left and right
    % 25 pieces of data per song
    for kk = 40:5:160
        test = song(1, Fs1*kk : Fs1*(kk+5));
%        test = resample(test, 20000, Fs1);
        vector = abs(spectrogram(test));
        vector = reshape(vector, [1, 8*32769]);
        A1 = [A1;vector];
    end
end
A1 = A1';
%%
[u_t1, s_t1, v_t1] = svd(A1 - mean(A1(:)), 'econ');
plot(diag(s_t1) ./ sum(diag(s_t1)), 'ro');
%%
vvv = v_t1';
plot3(vvv(2, 1:50), vvv(3, 1:50), vvv(4, 1:50), 'ro'); hold on;
plot3(vvv(2, 51:100), vvv(3, 51:100), vvv(4, 51:100), 'co'); 
plot3(vvv(2, 101:150), vvv(3, 101:150), vvv(4, 101:150), 'mo'); hold off;
legend('Classic', 'Rap', 'Rock')
xlabel('Mode 2'); ylabel('Mode 3'); zlabel('Mode 4')
%%
acct1_knn = []; acct1_nb = []; acct1_lda = [];
for test_trail = 1:1
    true = [ones(20,1); 2*ones(20,1); 3*ones(20,1)];
    q1 = randperm(50); q2 = randperm(50); q3 = randperm(50);
    xclas = v_t1(1:50, 2:4);
    xpop = v_t1(51:100, 2:4);
    xrock = v_t1(101:150, 2:4);
    xtrain_t1 = [xclas(q1(1:30), :); xpop(q2(1:30), :); xrock(q3(1:30),:)];
    xtest_t1 = [xclas(q1(31:end), :); xpop(q2(31:end), :); xrock(q3(31:end),:)];
    % knn
    ind = knnsearch(xtrain_t1, xtest_t1); 
    for i = 1:length(ind)
       if ind(i) <= 30
           ind(i) =  1;
       elseif ind(i) <= 60
           ind(i) = 2;
       else
           ind(i) = 3;
       end
    end
    temp = [ind==true];
    acct1_knn(test_trail) = sum(temp) / length(temp);
    subplot(3,3,1)
    bar(ind);
    title('kNN');
    xlabel('Test data'); ylabel('Label');
    
    % naive bayes
    ctrain = [ones(30,1); 2*ones(30,1); 3*ones(30,1)];
    nb = fitcnb(xtrain_t1, ctrain);
    pre = nb.predict(xtest_t1);
    temp = [pre== true];
    acct1_nb(test_trail) = sum(temp) / length(temp);
    subplot(3,3,2);
    bar(pre)
    title('Naive Bayes');
    xlabel('Test data'); ylabel('Label');
    
    % classify (Built in)
    pre = classify(xtest_t1, xtrain_t1, ctrain);
    temp = [pre== true];
    acct1_lda(test_trail) = sum(temp) / length(temp);
    subplot(3,3,3);
    bar(pre);
    title('LDA');
    xlabel('Test data'); ylabel('Label');    
end
%%
acct1_knn = mean(acct1_knn);
acct1_nb = mean(acct1_nb);
acct1_lda = mean(acct1_lda);
result = [acct1_knn;acct1_nb;acct1_lda];