% release resources and close figures
close all;
clear all;

% load files
adamRest=hhmbinread('hhm_konz1.hhm');
zsofiaRest=hhmbinread('hhmkonz2.hhm');
brigittaRest=hhmbinread('hhm_konz3.hhm');

% this transformation applies to the „2” signed HHMD 
% .ecg1 is a property of the object (for instance: adamRest)
adamRestMv = 3.3 / 4096 *(adamRest.ecg1 - 2048);
zsofiaRestMv = 3.3 / 4096 *(zsofiaRest.ecg1 - 2048);
brigittaRestMv = 3.3 / 4096 *(brigittaRest.ecg1 - 2048);

% this transformation applies to the „1” signed HHMD 
% adamRestMv = 3.3 / 8192 * (adamRest - 2048);
% zsofiaRestMv = 3.3 / 8192 *(zsofiaRest - 2048);
% brigittaRestMv = 3.3 / 8192 *(brigittaRest - 2048); 

% 100s
fs=100;

% visualize data
figure()
t=(1:length(adamRestMv))/fs; % time axis [sec]
plot(t,adamRestMv)

% visualize filtered data
% using 5th-order Butterworth filter
adamRestMvFiltered = butterworthFilter(adamRestMv,5);
figure()
t=(1:length(adamRestMvFiltered))/fs; % time axis [sec]
plot(t,adamRestMvFiltered)

function highPassFilteredSignal=butterworthFilter(ecg, order)
% ecg: hhm object's property
% order: Butterworth filter parameter, for instance: 5 -> 5th-order
% Butterworth filter
fs=1024; % sampled at this frequency [Hz]
n = order; % xth-order Butterworth filter parameter
fUpperCutOff = 0.5; % upper cutoff frequency [Hz]
fLowerCutOff = 20; % lower cutoff frequency [Hz]
% create two Butterwort filters, one for upper cutoff, one for lower cutoff
[bLower,aLower]=butter(n,fLowerCutOff/(fs*0.5),'low');
[bUpper,aUpper]=butter(n,fUpperCutOff/(fs*0.5),'high'); 
% filter signal
lowPassFilteredSignal = filtfilt(bLower, aLower, ecg);
highPassFilteredSignal = filtfilt(bUpper, aUpper, lowPassFilteredSignal);
end

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
