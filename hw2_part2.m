%% piano
clear all; close all; clc;
tr_piano=16; % record time in seconds
y1=audioread('music1.wav'); Fs1=length(y1)/tr_piano;
y1 = y1'/ 2;
%p8 = audioplayer(y1',Fs1); playblocking(p8);
L1 = tr_piano;
n1 = length(y1);
t1 = (1:length(y1))/Fs1;
k1=(2*pi/L1)*[0:n1/2-1 -n1/2:-1]; ks1=fftshift(k1);
freq1 = [];
Sgt_spec=[];
tslide=0:0.25:L1;
a = 55;
for j=1:length(tslide)
    g=exp(-a*(t1-tslide(j)).^2);
    Sg=g.*y1; 
    Sgt=fft(Sg);
    [V, I] = max(abs(Sgt));
    freq1 = [freq1; abs(k1(I))/(2*pi)];
    Sgt_spec=[Sgt_spec; abs(fftshift(Sgt)) / max(abs(Sgt))];
end

figure(1)
plot(tslide,freq1)
title('Center Frequency for the piano');
xlabel('Time(sec)');
ylabel('Frequency(Hertz)');

figure(2)
pcolor(tslide, ks1, abs(Sgt_spec).'), shading interp
set(gca,'Ylim',[1000 30000],'Fontsize',[14])
title('Spectrogram for the piano')
xlabel('Time(sec)');
ylabel('Frequency(\omega)');
colormap(hot)

%% record
tr_rec=14; % record time in seconds
y2=audioread('music2.wav'); Fs2=length(y2)/tr_rec;
y2 = y2';
%p8 = audioplayer(y2,Fs2); playblocking(p8);
n2 = length(y2);
t2 = (1:length(y2))/Fs2;
k2=(2*pi/(tr_rec))*[0:(n2/2-1) -n2/2:-1]; ks2=fftshift(k2);
freq2 = [];
Sgt_spec2=[];
tslide2=0:0.1:14;
a = 40;
for j=1:length(tslide2)
    g2=exp(-a*(t2-tslide2(j)).^2); % Gabor
    Sg2=g2.*y2; Sgt2=fft(Sg2);
    [V, I] = max(abs(Sgt2));
    freq2 = [freq2; abs(k2(I)/(2*pi))];
end

figure(3)
plot(tslide2, freq2);
title('Center Frequency for the recorder');
xlabel('Time(sec)');
ylabel('Frequency(Hertz)');

figure(4)
pcolor(tslide2, ks2, abs(Sgt_spec2).'), shading interp
title('Spectrogram for the recorder')
set(gca,'Ylim',[4000 30000],'Fontsize',[14])
xlabel('Time(sec)');
ylabel('Frequency(\omega)');
colormap(hot)