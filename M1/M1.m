% release resources and close figures
close all;
clear all;

%% Load files and set variables
adamRest=hhmbinread('rest/adam_relaxed.hhm');
zsofiaRest=hhmbinread('rest/zsofia_relaxed.hhm');
brigittaRest=hhmbinread('rest/brigitta_relaxed.hhm');

% % this transformation applies to the „2” signed HHMD 
% % .ecg1 is a property of the object (for instance: adamRest) and holds the
% % value of Einthoven I EKG signals
% adamRestMv = 3.3 / 4096 *(adamRest.ecg1 - 2048);
% zsofiaRestMv = 3.3 / 4096 *(zsofiaRest.ecg1 - 2048);
% brigittaRestMv = 3.3 / 4096 *(brigittaRest.ecg1 - 2048);
% 
% % this transformation applies to the „1” signed HHMD 
adamRestMv = 3.3 / 8192 * (adamRest.ecg1 - 2048);
zsofiaRestMv = 3.3 / 8192 *(zsofiaRest.ecg1 - 2048);
brigittaRestMv = 3.3 / 8192 *(brigittaRest.ecg1 - 2048); 
% 
fs=1000; % sampled at this frequency [Hz]
 
% % visualize data
% figure()
% t=(1:length(brigittaRestMv))/fs; % time axis [sec]
% plot(t,brigittaRestMv)

%% 1/a Filter and visualize signals

% filter signals using 5th-order Butterworth filter, 0.5 Hz lower and 40Hz
% upper cutoff frequencies
%adamRestMvFiltered = butterworthFilter(adamRestMv,5,40,0.5);
%zsofiaRestMvFiltered = butterworthFilter(zsofiaRestMv,5,40,0.5);
%brigittaRestMvFiltered = butterworthFilter(brigittaRestMv,5,40,0.5);

% visualize filtered signal
% figure()
% t=(1:length(brigittaRestMvFiltered))/fs; % time axis [sec]
% plot(t,brigittaRestMvFiltered)

% filter signals using 5th-order Butterworth filter, 0.5 Hz upper and 20Hz
% lower cutoff frequencies
adamRestMvFiltered = butterworthFilter(adamRestMv,5,20,0.5);
zsofiaRestMvFiltered = butterworthFilter(zsofiaRestMv,5,30,0.5);
brigittaRestMvFiltered = butterworthFilter(brigittaRestMv,5,20,0.5);

% % visualize filtered signal
% figure()
% t=(1:length(adamRestMvFiltered))/fs; % time axis [sec]
% subplot(3,1,1)
% plot(t,adamRestMv)
% title('Adam rest')
% subplot(3,1,2)
% plot(t,adamRestMvFiltered)
% 
% figure()
% t=(1:length(zsofiaRestMvFiltered))/fs; % time axis [sec]
% subplot(4,1,1)
% plot(t,zsofiaRestMv)
% title('Zsofia rest')
% subplot(4,1,2)
% plot(t,zsofiaRestMvFiltered)
% 
% figure()
% t=(1:length(brigittaRestMvFiltered))/fs; % time axis [sec]
% subplot(7,1,1)
% plot(t,brigittaRestMv)
% title('Brigitta rest')
% subplot(3,1,2)
% plot(t,brigittaRestMvFiltered)

% split signals
adamRestMvFilteredSplit = adamRestMvFiltered(50000:90000);
zsofiaRestMvFilteredSplit = zsofiaRestMvFiltered(40000:80000);
brigittaRestMvFilteredSplit = brigittaRestMvFiltered(10000:50000);

% visualize split signals
% t=(1:length(adamRestMvFilteredSplit))/fs; % all three signals have the same length by now
% figure()
% subplot(3,1,1) %(3 rows, 1th element)
% plot(t,adamRestMvFilteredSplit); title('Adam rest');xlabel('t [sec]');ylabel('U [mV]');
% subplot(3,1,2)
% plot(t,zsofiaRestMvFilteredSplit); title('Zsofia rest');xlabel('t [sec]');ylabel('U [mV]');
% subplot(3,1,3)
% plot(t,brigittaRestMvFilteredSplit); title('Brigitta rest');xlabel('t [sec]');ylabel('U [mV]');

%% 1/b Visualize in frequency domain

%frequencySpectrum(adamRestMvFilteredSplit,'Adam rest fft');
%frequencySpectrum(zsofiaRestMvFilteredSplit,'Zsofia rest fft');
%frequencySpectrum(brigittaRestMvFilteredSplit,'Brigitta rest fft');
% 
% %holding breath for 30 seconds: Adam 
% adamHoldBreath30Sec=hhmbinread('hold_breath/adam_hold_breath_30sec.hhm');
% adamHoldBreath30SecFiltered = butterworthFilter(adamHoldBreath30Sec.ecg1 - mean(adamHoldBreath30Sec.ecg1),5,20,0.5);
% frequencySpectrum(adamHoldBreath30SecFiltered,'Adam - holding breath for 30 seconds - fft'); %visualize in frequency domain 
%  
% %holding breath for 30 seconds: Zsofia 
% zsofiaHoldBreath30Sec=hhmbinread('hold_breath/zsofia_hold_breath_30sec.hhm');
% zsofiaHoldBreath30SecFiltered = butterworthFilter(zsofiaHoldBreath30Sec.ecg1 - mean(zsofiaHoldBreath30Sec.ecg1),5,20,0.5);
% frequencySpectrum(zsofiaHoldBreath30SecFiltered,'Zsofia - holding breath for 30 seconds - fft'); %visualize in frequency domain 

% %holding breath for 30 seconds: Brigitta
 brigittaHoldBreath30Sec=hhmbinread('hold_breath/brigitta_hold_breath_30sec.hhm');
% brigittaHoldBreath30SecFiltered = butterworthFilter(brigittaHoldBreath30Sec.ecg1 - mean(brigittaHoldBreath30Sec.ecg1),5,20,0.5);
% frequencySpectrum(brigittaHoldBreath30SecFiltered,'Brigitta - holding breath for 30 seconds - fft'); %visualize in frequency domain 

% %controlling breathing: Adam 
% adamControlBreath=hhmbinread('control_breath/adam_control_breath.hhm');
% adamControlBreathFiltered = butterworthFilter(adamControlBreath.ecg1 - mean(adamControlBreath.ecg1),5,20,0.5);
% frequencySpectrum(adamControlBreathFiltered,'Adam - controlled breathing - fft'); %visualize in frequency domain 
%  
% %controlling breathing: Zsofia 
% zsofiaControlBreath=hhmbinread('control_breath/zsofia_control_breath.hhm');
% zsofiaControlBreathFiltered = butterworthFilter(zsofiaControlBreath.ecg1 - mean(zsofiaControlBreath.ecg1),5,20,0.5);
% frequencySpectrum(zsofiaControlBreathFiltered,'Zsofia - controlled breathing - fft'); %visualize in frequency domain 
% 
%controlling breathing: Brigitta
 brigittaControlBreath=hhmbinread('control_breath/brigitta_control_breath.hhm');
% brigittaControlBreathFiltered = butterworthFilter(brigittaControlBreath.ecg1 - mean(brigittaControlBreath.ecg1),5,20,0.5);
% frequencySpectrum(brigittaControlBreathFiltered,'Brigitta - controlled breathing - fft'); %visualize in frequency domain
% 
% %Physical activity: Brigitta
 brigittaPhysicalActivity=hhmbinread('squat/brigitta_squat.hhm');
% brigittaPhysicalActivityFiltered = butterworthFilter(brigittaPhysicalActivity.ecg1 - mean(brigittaPhysicalActivity.ecg1),5,20,0.5);
% frequencySpectrum(brigittaPhysicalActivityFiltered,'Brigitta - physical activity - fft'); %visualize in frequency domain

% %% 2 task

fs=1000; % sampled at this frequency [Hz]

% load red ppgl signal into variables
brigittaRestPpglRed=brigittaRest.ppgl_red;
brigittaRestPpglRedSplit = brigittaRestPpglRed(10000:length(brigittaRestPpglRed)); % be aware of matrix size

brigittaHoldBreathPpglRed = brigittaHoldBreath30Sec.ppgl_red;
brigittaHoldBreathPpglRedSplit = brigittaHoldBreathPpglRed(10000:length(brigittaHoldBreathPpglRed));

brigittaSquat = hhmbinread('squat/brigitta_squat.hhm'); % this file was not loaded
brigittaSquatPpglRed = brigittaSquat.ppgl_red;
brigittaSquatPpglRedSplit = brigittaSquatPpglRed(10000:length(brigittaSquatPpglRed));

brigittaControlBreathPpglRed = brigittaControlBreath.ppgl_red;
brigittaControlBreathPpglRedSplit = brigittaControlBreathPpglRed(10000:length(brigittaControlBreathPpglRed));

% visualize signals
figure()

t=(1:length(brigittaRestPpglRedSplit))/fs;
subplot(3,1,1); 
plot(t,3.3 / 8192 *(brigittaRestPpglRedSplit - 2048)); title('brigitta rest ppgl red'); % same rescaling used as in the beginning, be aware of HHMD version
xlabel('t [sec]');ylabel('U [mV]');

t=(1:length(brigittaHoldBreathPpglRedSplit))/fs;
subplot(3,1,2); 
plot(t,3.3 / 8192 *(brigittaHoldBreathPpglRedSplit - 2048)); title('brigitta holding breath ppgl red');
xlabel('t [sec]');ylabel('U [mV]');

figure()

t=(1:length(brigittaControlBreathPpglRedSplit))/fs;
subplot(4,1,1); 
plot(t,3.3 / 8192 *(brigittaControlBreathPpglRedSplit - 2048)); title('brigitta controlling breathing ppgl red');
xlabel('t [sec]');ylabel('U [mV]');

t=(1:length(brigittaSquatPpglRedSplit))/fs;
subplot(4,1,2); 
plot(t,3.3 / 8192 *(brigittaSquatPpglRedSplit - 2048)); title('brigitta squat ppgl red');
xlabel('t [sec]');ylabel('U [mV]');

% % load infrared ppgl signal into variables 
% adamRestPpglInfrared=adamRest.ppgl_nir;
% adamRestPpglInfraredSplit = adamRestPpglInfrared(1:30000); % be aware of matrix size
% 
% adamHoldBreathPpglInfrared = adamHoldBreath30Sec.ppgl_nir;
% adamHoldBreathPpglInfraredSplit = adamHoldBreathPpglInfrared(1:30000);
% 
% adamControlBreathPpglInfrared = adamControlBreath.ppgl_nir;
% adamControlBreathPpglInfraredSplit = adamControlBreathPpglInfrared(1:30000);
% 
% adamSquat = hhmbinread('squat/adam_squat.hhm'); % this file was not loaded
% adamSquatPpglInfrared = adamSquat.ppgl_nir;
% adamSquatPpglInfraredSplit = adamSquatPpglInfrared(1:30000);

% % visualize signals
% figure()
% 
% t=(1:length(brigittaRestPpglInfraredSplit))/fs;
% subplot(5,1,1); 
% plot(t,3.3 / 4096 *(brigittaRestPpglInfraredSplit - 2048)); title('brigitta rest ppgl infrared'); % same rescaling used as in the beginning, be aware of HHMD version
% xlabel('t [sec]');ylabel('U [mV]');
% 
% t=(1:length(brigittaHoldBreathPpglInfraredSplit))/fs;
% subplot(5,1,2); 
% plot(t,3.3 / 4096 *(brigittaHoldBreathPpglInfraredSplit - 2048)); title('brigitta holding breath ppgl infrared');
% xlabel('t [sec]');ylabel('U [mV]');
% 
% t=(1:length(brigittaControlBreathPpglInfraredSplit))/fs;
% subplot(6,1,1); 
% plot(t,3.3 / 4096 *(brigittaControlBreathPpglInfraredSplit - 2048)); title('brigitta controlling breathing ppgl infrared');
% xlabel('t [sec]');ylabel('U [mV]');
% 
% t=(1:length(brigittaSquatPpglInfraredSplit))/fs;
% subplot(6,1,2); 
% plot(t,3.3 / 4096 *(brigittaSquatPpglInfraredSplit - 2048)); title('brigitta squat ppgl infrared');
% xlabel('t [sec]');ylabel('U [mV]');

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
% 
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
