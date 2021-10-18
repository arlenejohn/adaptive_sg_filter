%% Adaptive SG filtering Example Program


clear all
close all
clc

%%parameter definition
Mmax=20; % Maximum window size for window varying denoising
pmax=5; % Maximum order for window varying denoising
M=15; % general window size used
p=3; % general order used


%%%% loading the signal (here a clean signal is loaded and noise is added
%%%% to it. However, in real-life scenarios, the signal is just loaded here

load('aami3am_noisy_L.mat')
y1=noisy;

%%% normalize the signal
y1=y1-mean(y1);
sig=y1/max(abs(y1)); %normalizing the signal
noisy=sig;

%OR JUST LOAD THE NOISY SIGNAL

     
%% SG-filter based Denoising
type='L'; %specify noist type as gaussian, laplacian, or uniform for sigma estimate calculation

% here, the G-O-R algorithm is executed. You could choose between G-FL, G-FL-R,
% G-O, OR G-O-R by uncommenting the corresponding lines. G-O-R is the best performing algorithm. Hence is used
% here.

%[denoised_signal,window_order] = den_win(Mmax,p,noisy,type); %window varying denoising (G-FL)
%[denoised_signal,window_order] = den_win_reg(Mmax,p,noisy,type);%regularized window varying denoising(G-FL-R)
%[denoised_signal,window_order] = den_ord(M,pmax,noisy,type); %order varying denoising (G-O)
[denoised_signal,window_order] = den_ord_reg(M,pmax,noisy,type);%regularized order varying denoising(G-O-R)


%%%%% plotting when the clean signal is not available
figure
subplot(3,1,1)
plot(noisy,'k')
legend('noisy/ input')
grid on
subplot(3,1,2)
plot(noisy,'k')
hold on
plot(denoised_signal,'r')
legend('noisy/ input','denoised')
grid on
subplot(3,1,3)
plot(window_order,'m')
legend('order/window')
grid on

