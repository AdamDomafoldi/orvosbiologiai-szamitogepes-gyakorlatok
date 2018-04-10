% release resources and close figures
close all;
clear all;

%% 3 task

adamBloodPressure=hhmbinread('blood_pressure/adam.hhm');
zsofiaBloodPressure=hhmbinread('blood_pressure/zsofia.hhm');
brigittaBloodPressure=hhmbinread('blood_pressure/brigitta.hhm');

% set timescale

fs=1000; 

tAdam=(1:length(adamBloodPressure))/fs;
tZsofia=(1:length(zsofiaBloodPressure))/fs;
tBrigitta=(1:length(brigittaBloodPressure))/fs;

% visualize blood pressure and PPG connection
figure();

subplot(3,1,1); 
plot(tAdam,((adamBloodPressure.press-175)/15.5),'red');
title('Adam: blood pressure and PPG connection');
hold on
plot(tAdam,200*(adamBloodPressure.press-175)/15.5,'black'); % zooming
xlabel('t [sec]');

subplot(3,1,2); 
plot(tZsofia,((zsofiaBloodPressure.press-175)/15.5),'red');
title('Zsofia: blood pressure and PPG connection');
hold on
plot(tZsofia,200*(zsofiaBloodPressure.press-175)/15.5,'black'); % zooming
xlabel('t [sec]');
ylabel('blood pressure [mmHg]');

subplot(3,1,3); 
plot(tBrigitta,((brigittaBloodPressure.press-175)/15.5),'red');
title('Brigitta: blood pressure and PPG connection');
hold on
plot(tBrigitta,200*(brigittaBloodPressure.press-175)/15.5,'black'); % zooming
xlabel('t [sec]');

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