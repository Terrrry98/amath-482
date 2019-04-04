clear all; close all; clc;
load handel
v = y'/2;
v = v(1:73112);
L = 9;
n = length(v);
t = (1:length(v))/Fs;
k=(2*pi/(L))*[0:(n-1)/2 -(n)/2:-1]; ks=fftshift(k);

%% 1.Gaussian. width = 1, delta_t = 0.1
vgt_spec=[];
tslide=0:0.1:9;
a = 1;
for j=1:length(tslide)
    g=exp(-a*(t-tslide(j)).^2); 
    vg=g.*v; vgt=fft(vg);
    vgt_spec=[vgt_spec; abs(fftshift(vgt))];
end
pcolor(tslide, ks, vgt_spec.'), shading interp
xlabel('Time(sec)');
ylabel('Frequency(\omega)');
title('Spectrogram through Gabor filtering with width = 1');
colormap(hot)


%% 2. Gaussian. width = 50, delta_t = 0.1
close all;
vgt_spec=[];
tslide=0:0.1:9;
a = 50;
for j=1:length(tslide)
    g=exp(-a*(t-tslide(j)).^2); 
    vg=g.*v; vgt=fft(vg);
    vgt_spec=[vgt_spec; abs(fftshift(vgt))];
end
pcolor(tslide,ks,vgt_spec.'), shading interp
xlabel('Time(sec)');
ylabel('Frequency(\omega)');
title('Spectrogram through Gabor filtering with width = 50');
colormap(hot)

%% 2. Gaussian. width = 0.01, delta_t = 0.1
close all;
vgt_spec=[];
tslide=0:0.1:9;
a = 0.01;
for j=1:length(tslide)
    g=exp(-a*(t-tslide(j)).^2); 
    vg=g.*v; vgt=fft(vg);
    vgt_spec=[vgt_spec; abs(fftshift(vgt))];
end
pcolor(tslide,ks,vgt_spec.'), shading interp
xlabel('Time(sec)');
ylabel('Frequency(\omega)');
title('Spectrogram through Gabor filtering with width = 0.01');
colormap(hot)

%% 2. Gaussian. width = 10, delta_t = 0.1
close all;
vgt_spec=[];
tslide=0:0.1:9;
a = 10;
for j=1:length(tslide)
    g=exp(-a*(t-tslide(j)).^2); 
    vg=g.*v; vgt=fft(vg);
    vgt_spec=[vgt_spec; abs(fftshift(vgt))];
end
pcolor(tslide,ks,vgt_spec.'), shading interp
xlabel('Time(sec)');
ylabel('Frequency(\omega)');
title('Spectrogram through Gabor filtering with width = 10');
colormap(hot)


%% 3. Gaussian. width = 1, delta_t = 0.05
close all;
vgt_spec=[];
tslide=0:0.05:9;
a = 1;
for j=1:length(tslide)
    g=exp(-a*(t-tslide(j)).^2);
    vg=g.*v; vgt=fft(vg);
    vgt_spec=[vgt_spec; abs(fftshift(vgt))];
end
pcolor(tslide,ks,vgt_spec.'), shading interp
xlabel('Time(sec)');
ylabel('Frequency(\omega)');
title('Spectrogram through Gabor filtering with time step = 0.05');
colormap(hot)

%% 3. Gaussian. width = 1, delta_t = 0.5
close all;
vgt_spec=[];
tslide=0:0.5:9;
a = 1;
for j=1:length(tslide)
    g=exp(-a*(t-tslide(j)).^2);
    vg=g.*v; vgt=fft(vg);
    vgt_spec=[vgt_spec; abs(fftshift(vgt))];
end
pcolor(tslide,ks,vgt_spec.'), shading interp
xlabel('Time(sec)');
ylabel('Frequency(\omega)');
title('Spectrogram through Gabor filtering with time step = 0.5');
colormap(hot)

%% 3. Gaussian. width = 1, delta_t = 1
close all;
vgt_spec=[];
tslide=0:1:9;
a = 1;
for j=1:length(tslide)
    g=exp(-a*(t-tslide(j)).^2);
    vg=g.*v; vgt=fft(vg);
    vgt_spec=[vgt_spec; abs(fftshift(vgt))];
end
figure(2)
pcolor(tslide,ks,vgt_spec.'), shading interp
xlabel('Time(sec)');
ylabel('Frequency(\omega)');
title('Spectrogram through Gabor filtering with time step = 1');
colormap(hot)

%% 4:
% Gaussian filter
close all;
vgt_spec=[];
tslide=0:0.1:9;
a  = 5;
for j=1:length(tslide)
    g= exp(-a.*(t - tslide(j)).^2);
    vg=g.*v;
    vgt = fft(vg);
    vgt_spec=[vgt_spec; abs(fftshift(vgt))];
end
subplot(3,1,1);
pcolor(tslide,ks,vgt_spec.'), shading interp
xlabel('Time(sec)');
ylabel('Frequency(\omega)');
title('Gaussian window')
colormap(hot)

% step function
vgt_spec=[];
tslide=0:0.1:9;
width = 1;

for j=1:length(tslide)
    s = (abs(t - tslide(j)) < width);
    vg=s.*v;
    vgt = fft(vg);
    vgt_spec=[vgt_spec; abs(fftshift(vgt))];
    %subplot(3,1,1), plot(t,v,'k',t,s,'r')
    %subplot(3,1,2), plot(t,vg,'k')
    %subplot(3,1,3), plot(ks,abs(fftshift(vgt))/max(abs(vgt)))
    %pause(0.1)
end
subplot(3,1,2);
pcolor(tslide,ks,vgt_spec.'), shading interp
xlabel('Time(sec)');
ylabel('Frequency(\omega)');
title('Shannon window');
colormap(hot)

% Mexican hat wavelet
vgt_spec=[];
tslide=0:.1:9;
a = 0.5;
for j=1:length(tslide)
    m = 2/(sqrt(3*a)*(pi)^(1/4))*(1-((t - tslide(j))/a).^2)...
        .* exp(-(t - tslide(j)).^2 / (2*a^2));
    vg=m.*v; vgt = fft(vg);
    vgt_spec=[vgt_spec; abs(fftshift(vgt))];
end
subplot(3,1,3)
pcolor(tslide,ks,vgt_spec.'), shading interp
xlabel('Time(sec)');
ylabel('Frequency(\omega)');
title('Mexican hat wavelet');
colormap(hot)


%% plot 3 filters in the same figure
j = 45;
g= exp(-5.*(t - tslide(j)).^2);
s = (abs(t - tslide(j)) < width);
m = 2/(sqrt(3*a)*(pi)^(1/4))*(1-((t - tslide(j))/a).^2)...
        .* exp(-(t - tslide(j)).^2 / (2*a^2));
figure(2);

subplot(3,1,1), plot(t,v,'k',t,g,'r');
title('Gaussian filter with a = 5');
xlabel('Time [sec]');
ylabel('Amplitude');
subplot(3,1,2), plot(t,v,'k',t,s,'r');
xlabel('Time [sec]');
ylabel('Amplitude');
title('Shannon window with width = 1');
subplot(3,1,3), plot(t,v,'k',t,m,'r');
xlabel('Time [sec]');
ylabel('Amplitude');
title('Mexican hat wavelet with \sigma = 0.5');

