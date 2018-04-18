% release resources and close figures
close all;
clear all;

%% Load files and set variables

signal=hhmbinread('mm_nyugalmi.hhm');

ecg = toMillivolt(signal.ecg1);

fs=1000;
ecg=ecg(25*fs:end); %kezdo tranziens levagasa
t=(1:length(ecg))/fs;

[baseLineNum, baseLineDum]=butter(3,0.5/500);
figure;
freqz(baseLineNum, baseLineDum);

figure();
plot(t,ecg);
figure;
plot(t,ecg);

%50/500=0.1 az 50Hz
[f50num, f50den]=butter(3,[0.09 0.11],'stop');
% [f50num f50den]=cheby1(3,3,[0.05 0.15],'stop'); %pass band ripple 3dB 
figure();
freqz(f50num, f50den); %frekvencia valasz abrazolasa
ecg50=filtfilt(f50num,f50den,ecg);
figure();
plot(t,ecg50);

% 35%500=0.07
[lp35num, lp35den]=butter(3,0.07);
figure();
freqz(lp35num, lp35den); %frekvencia valasz abrazolasa
ecglp35=filter(lp35num, lp35den,ecg50);
figure();
plot(t,ecglp35);

%qrsDetect(ecg);
%normalo(ecg50);
% 
function out = toMillivolt(ECGsignal)
    out = 3.3 / 8192 * (ECGsignal - 2048);
end





