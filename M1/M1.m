% release resources and close figures
close all;
clear all;

%% Load files and set variables
adamRest=hhmbinread('rest/adam_relaxed.hhm');
zsofiaRest=hhmbinread('rest/zsofia_relaxed.hhm');
brigittaRest=hhmbinread('rest/brigitta_relaxed.hhm');

% this transformation applies to the „1” signed HHMD 
% .ecg1 is a property of the object (for instance: adamRest) and holds the
% value of Einthoven I EKG signals
adamRestMv = toMillivolt(adamRest.ecg1);
zsofiaRestMv = toMillivolt(zsofiaRest.ecg1);
brigittaRestMv = toMillivolt(brigittaRest.ecg1); 

fs=1000; % sampled at this frequency [Hz]
 
% % visualize data
% % this plots the raw signal
% figure()
% t=(1:length(brigittaRestMv))/fs; % time axis [sec]
% plot(t,brigittaRestMv)

%% 1/a Filter and visualize signals

% filter signals using 5th-order Butterworth filter, 0.5 Hz upper and
% 20Hz-30Hz lower cutoff frequencies
adamRestMvFiltered = butterworthFilter(adamRestMv,5,20,0.5);
zsofiaRestMvFiltered = butterworthFilter(zsofiaRestMv,5,30,0.5);
brigittaRestMvFiltered = butterworthFilter(brigittaRestMv,5,20,0.5);

% % visualize filtered signal
% figure()
% t=(1:length(adamRestMvFiltered))/fs; % time axis [sec]
% subplot(2,1,1)
% plot(t,adamRestMv)
% title('Adam rest - filtered EKG signal')
% subplot(2,1,2)
% plot(t,adamRestMvFiltered)
% 
% figure()
% t=(1:length(zsofiaRestMvFiltered))/fs; % time axis [sec]
% subplot(3,1,1)
% plot(t,zsofiaRestMv)
% title('Zsofia rest - filtered EKG signal')
% subplot(3,1,2)
% plot(t,zsofiaRestMvFiltered)
% 
% figure()
% t=(1:length(brigittaRestMvFiltered))/fs; % time axis [sec]
% subplot(4,1,1)
% plot(t,brigittaRestMv)
% title('Brigitta rest - filtered EKG signal')
% subplot(4,1,2)
% plot(t,brigittaRestMvFiltered)

% split signals
adamRestMvFilteredSplit = adamRestMvFiltered(50000:90000);
zsofiaRestMvFilteredSplit = zsofiaRestMvFiltered(40000:80000);
brigittaRestMvFilteredSplit = brigittaRestMvFiltered(10000:50000);

% visualize split signals
t=(1:length(adamRestMvFilteredSplit))/fs; % all three signals have the same length (40 sec) by now

figure()
subplot(5,1,1) %(3 rows, 1th element)
plot(t,adamRestMvFilteredSplit);
title('Adam rest - filtered and split EKG signal');xlabel('t [sec]');ylabel('U [mV]');

subplot(5,1,2)
plot(t,zsofiaRestMvFilteredSplit);
title('Zsofia rest - filtered and split EKG signal');xlabel('t [sec]');ylabel('U [mV]');

subplot(5,1,3)
plot(t,brigittaRestMvFilteredSplit);
title('Brigitta rest - filtered and split EKG signal');xlabel('t [sec]');ylabel('U [mV]');

%% 1/b Visualize in frequency domain
% fast fourier transform ekg signals and visualize them, original unfiltered
% signal is used, no split
frequencySpectrum(adamRestMv,'Adam rest - EKG frequency spectrum (original length)');
frequencySpectrum(zsofiaRestMv,'Zsofia rest - EKG frequency spectrum (original length)');
frequencySpectrum(brigittaRestMv,'Brigitta rest - EKG frequency spectrum (original length)'); 

%% 2/a
% ekg signals
% holding breath for 30 seconds: Brigitta
brigittaHoldBreath30Sec=hhmbinread('hold_breath/brigitta_hold_breath_30sec.hhm');
brigittaHoldBreath30SecFiltered = butterworthFilter(brigittaHoldBreath30Sec.ecg1,5,20,0.5);
frequencySpectrum(brigittaHoldBreath30SecFiltered,'Brigitta - holding breath for 30 seconds - filtered EKG frequency spectrum (original length)'); %visualize in frequency domain 

% controlling breathing: Brigitta
brigittaControlBreath=hhmbinread('control_breath/brigitta_control_breath.hhm');
brigittaControlBreathFiltered = butterworthFilter(brigittaControlBreath.ecg1,5,20,0.5);
frequencySpectrum(brigittaControlBreathFiltered,'Brigitta - controlled breathing - filtered EKG frequency spectrum (original length)'); %visualize in frequency domain

% physical activity: Brigitta
brigittaPhysicalActivity=hhmbinread('squat/brigitta_squat.hhm');
brigittaPhysicalActivityFiltered = butterworthFilter(brigittaPhysicalActivity.ecg1,5,20,0.5);
frequencySpectrum(brigittaPhysicalActivityFiltered,'Brigitta - physical activity - filtered EKG frequency spectrum (original length)'); %visualize in frequency domain

% ppgl signals
fs=1000; % sampled at this frequency [Hz]

%load red ppgl RED signals into variables and split them
brigittaRestPpglRed=brigittaRest.ppgl_red;
brigittaRestPpglRedSplit = brigittaRestPpglRed(10000:length(brigittaRestPpglRed));

brigittaHoldBreathPpglRed = brigittaHoldBreath30Sec.ppgl_red;
brigittaHoldBreathPpglRedSplit = brigittaHoldBreathPpglRed(10000:length(brigittaHoldBreathPpglRed));

brigittaPhysicalActivityPpglRed = brigittaPhysicalActivity.ppgl_red;
brigittaPhysicalActivityPpglRedSplit = brigittaPhysicalActivityPpglRed(10000:length(brigittaPhysicalActivityPpglRed));

brigittaControlBreathPpglRed = brigittaControlBreath.ppgl_red;
brigittaControlBreathPpglRedSplit = brigittaControlBreathPpglRed(10000:length(brigittaControlBreathPpglRed));

% visualize ppgl RED signals
figure();

t=(1:length(brigittaRestPpglRedSplit))/fs;
subplot(6,1,1); 
plot(t,toMillivolt(brigittaRestPpglRedSplit));
title('Brigitta rest - red split ppgl signal');
xlabel('t [sec]');
ylabel('U [mV]');

t=(1:length(brigittaHoldBreathPpglRedSplit))/fs;
subplot(6,1,2); 
plot(t,toMillivolt(brigittaHoldBreathPpglRedSplit));
title('Brigitta holding breath - red split ppgl signal');
xlabel('t [sec]');ylabel('U [mV]');

figure();

t=(1:length(brigittaControlBreathPpglRedSplit))/fs;
subplot(7,1,1); 
plot(t,toMillivolt(brigittaControlBreathPpglRedSplit));
title('Brigitta controlling breathing - red split ppgl signal');
xlabel('t [sec]');
ylabel('U [mV]');

t=(1:length(brigittaPhysicalActivityPpglRedSplit))/fs;
subplot(7,1,2); 
plot(t,toMillivolt(brigittaPhysicalActivityPpglRedSplit));
title('Brigitta squat - red split ppgl signal');
xlabel('t [sec]');
ylabel('U [mV]');

% load INFRARED (NIR) ppgl signals into variables 
brigittaRestPpglNir=brigittaRest.ppgl_nir;
brigittaRestPpglNirSplit = brigittaRestPpglNir(10000:length(brigittaRestPpglNir));

brigittaHoldBreathPpglNir = brigittaHoldBreath30Sec.ppgl_nir;
brigittaHoldBreathPpglNirSplit = brigittaHoldBreathPpglNir(10000:length(brigittaHoldBreathPpglNir));
 
brigittaPhysicalActivityPpglNir = brigittaPhysicalActivity.ppgl_nir;
brigittaPhysicalActivityPpglNirSplit = brigittaPhysicalActivityPpglNir(10000:length(brigittaPhysicalActivityPpglNir));
 
brigittaControlBreathPpglNir = brigittaControlBreath.ppgl_nir;
brigittaControlBreathPpglNirSplit = brigittaControlBreathPpglNir(10000:length(brigittaControlBreathPpglNir));

% visualize ppgl INFRARED (NIR) signals
figure();

t=(1:length(brigittaRestPpglNirSplit))/fs;
subplot(6,1,1); 
plot(t,toMillivolt(brigittaRestPpglRedSplit));
title('Brigitta rest - Nir split ppgl signal');
xlabel('t [sec]');
ylabel('U [mV]');

t=(1:length(brigittaHoldBreathPpglNirSplit))/fs;
subplot(6,1,2); 
plot(t,toMillivolt(brigittaHoldBreathPpglNirSplit));
title('Brigitta holding breath - Nir split ppgl signal');
xlabel('t [sec]');ylabel('U [mV]');

figure();

t=(1:length(brigittaControlBreathPpglNirSplit))/fs;
subplot(7,1,1); 
plot(t,toMillivolt(brigittaControlBreathPpglNirSplit));
title('Brigitta controlling breathing - Nir split ppgl signal');
xlabel('t [sec]');
ylabel('U [mV]');

t=(1:length(brigittaPhysicalActivityPpglNirSplit))/fs;
subplot(7,1,2); 
plot(t,toMillivolt(brigittaPhysicalActivityPpglNirSplit));
title('Brigitta squat - Nir split ppgl signal');
xlabel('t [sec]');
ylabel('U [mV]');

%% 2/b

% fast fourier transform ppgl RED signals and visualize them, original unfiltered signal is used, no split
frequencySpectrum(brigittaControlBreathPpglRed,'Brigitta control breathing - ppgl RED frequency spectrum (original length)');
frequencySpectrum(brigittaHoldBreathPpglRed,'Brigitta hold breath - ppgl RED frequency spectrum (original length)');
frequencySpectrum(brigittaRestPpglRed,'Brigitta rest - ppgl RED frequency spectrum (original length)'); 
frequencySpectrum(brigittaPhysicalActivityPpglRed,'Brigitta physical activity - ppgl RED frequency spectrum (original length)'); 

% fast fourier transform ppgl NIR signals and visualize them, original unfiltered signal is used, no split
frequencySpectrum(brigittaControlBreathPpglNir,'Brigitta control breathing - ppgl NIR frequency spectrum (original length)');
frequencySpectrum(brigittaHoldBreathPpglNir,'Brigitta hold breath - ppgl NIR frequency spectrum (original length)');
frequencySpectrum(brigittaRestPpglNir,'Brigitta rest - ppgl NIR frequency spectrum (original length)'); 
frequencySpectrum(brigittaPhysicalActivityPpglNir,'Brigitta physical activity - ppgl NIR frequency spectrum (original length)'); 

function out = toMillivolt(ECGsignal)
    out = 3.3 / 8192 * (ECGsignal - 2048);
end

%% Frequency spectrum
function frequencySpectrum(ecg,titleOfDiagram)
    % @ecg: hhm object's property
    % @titleOfDiagram: title of the diagram
    fs=1000; % sampled at this frequency [Hz]   
    Y = fft(ecg); % fast fourier transformation
    L=length(ecg); % signal length
    signal = abs(Y/L); % absolute value of frequency
    spektrum = signal(1:L/2+1); % cut the signal into half
    f = fs *(0:(L/2))/L; % frequency axsis [Hz] 
    % visualize half signal
    figure();
    plot(f,spektrum) 
    title(titleOfDiagram)
    xlabel('f [Hz]')
    ylabel('amplification')
    xlim([0 10]) % X axis
end

%% Butterworth filter
function highPassFilteredSignal=butterworthFilter(ecg, order, lowerCutOff,upperCutOff)
    % @ecg: hhm object's property
    % @order: Butterworth filter parameter, for instance: 5 -> 5th-order
    % Butterworth filter
    % @lowerCutOff: lower cutoff frequency [Hz]
    % @upperCutOff: upper cutoff frequency [Hz]
    fs=1024; % sampled at this frequency [Hz]
    n = order; % xth-order Butterworth filter parameter
    fUpperCutOff = upperCutOff; % upper cutoff frequency [Hz]
    fLowerCutOff = lowerCutOff; % lower cutoff frequency [Hz]
    % create two Butterwort filters, one for upper cutoff, one for lower cutoff
    [bLower,aLower]=butter(n,fLowerCutOff/(fs*0.5),'low');
    [bUpper,aUpper]=butter(n,fUpperCutOff/(fs*0.5),'high'); 
    % filter signal
    lowPassFilteredSignal = filtfilt(bLower, aLower, ecg);
    highPassFilteredSignal = filtfilt(bUpper, aUpper, lowPassFilteredSignal);
    end

%% HHM file reader
function o=hhmbinread(filename)
    %HHM file reader
    % Input: name of HHM file (with full path and extension)
    % Output: structure with fields:
    %  ecg1: ECG lead Einthoven I.
    %  ecg2: ECG lead Einthoven II.
    %  press: cuff pressure captured from analog/digital converter. Offset is
    %     not compensated. 15.5 LSB/mmHg
    %  ppgl_red: PPG left red
    %  ppgl_red_dc: PPG left red without DC level filtering (smaller gain)
    %  ppgr_nir: PPG right near infra red
    %  ppgl_nir_dc: PPG left nir without DC level filtering (smaller gain)
    %  ppgl_nir: PPG left near infra red
    %
    % Every channel is sampled with 1 kHz
    %Csordás Péter 2006
    fid=fopen(filename,'r','b');
    if fid==-1 error('File not found'); end

    buffer=fread(fid,'bit12=>double'); 
    idx=find(buffer<0);
    if(~isempty(idx))
        buffer(idx)=buffer(idx)+4096; %overflow correction
    end

    o.ecg1=buffer(1:8:end);
    o.ecg2=buffer(2:8:end);
    o.press=buffer(3:8:end);
    o.ppgl_red=buffer(4:8:end);
    o.ppgl_red_dc=buffer(5:8:end);
    o.ppgr_nir=buffer(6:8:end);
    o.ppgl_nir_dc=buffer(7:8:end);
    o.ppgl_nir=buffer(8:8:end);
end
