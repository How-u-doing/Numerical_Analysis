% Signal processing example of addtion of 2 sine waves
% According to Data4Bio, YouTube lectures by University of Washington

L=1;        % lasting time
dt=0.001;   % time interval
t=0:dt:L;   % time vector

x=sin(2*pi*50*t)+sin(2*pi*120*t); % pure source sound

figure(1);
subplot(2,1,1);
plot(t,x,'b','LineWidth',1.2);
hold on

% add noise
noise_coef=2.5;
y=x + noise_coef * randn(size(x));
plot(t,y,'r');
axis([0 .25 -5 5]);
legend('Clean','Noise');

% compute the (Fast) Fourier Transform
figure(2);
N=length(t);
Y=fft(y,N); % compute the FFT data points
PSD=Y.*conj(Y)/N;   % Power Spectrum Density...
% how much power in each frequency

freq=1/L*(0:N);
L1=1:floor(N/2);
plot(freq(L1),PSD(L1));
xlabel('Frequency (Hz)')
title('Power Spectrum Density')

% use PSD to filter out noise
% see the graph to adjust the threshold
idx=PSD>50; % find all indices with large power than threshold (it can be 50)  
PSD=PSD.*idx;
hold on
plot(freq(L1),PSD(L1));
legend('Original','Filtered')


Y=idx.*Y;
yfilt=ifft(Y); % inverse FFT & get filtered data

% go back to figure 1 and add filtered signal
figure(1)
subplot(2,1,2)
plot(t,x,'b','LineWidth',1.2);
hold on
plot(t,yfilt,'r');
axis([0 .25 -5 5]);
legend('Clean','Filtered');


