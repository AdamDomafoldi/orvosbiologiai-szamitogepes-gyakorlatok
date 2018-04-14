% release resources and close figures
close all;
clear all;

%% 3 task

% load files and set variables
signal=hhmbinread('blood_pressure/brigitta.hhm');
ecgMv = toMillivolt(signal.ecg1); % get ecg1 and covert it to millivolt
ppglNir = signal.ppgl_nir; % get nir value
pressure = toPressure(signal.press); % get pressure value and convert it

fs = 1000;
t = (0:length(ecgMv)-1)/fs;

% filter signals
% zero-phase digital filtering
[b,a] = butter(5,0.044);
ppglNirFilted = filtfilt(b,a,ppglNir);
% butterworth filtering
ecgMvFiltered = butterworthFilter(ecgMv,5,20,0.5);

% visualize data
figure();
plot(t,ppglNirFilted/20,'r')
title('PPG - Blood pressure');
xlabel('t [ms]');
ylabel(' Blood pressure [mmHg]');
hold on;
plot(t,pressure,'g');
plot(t,ecgMv*100,'b');
legend('PPG','Blood pressure','EKG');

hold off;

% PPG peak
[var, ppgBeginning] = findpeaks(ppglNirFilted * -1, 'minpeakdistance',600);
figure();
plot(t, ppglNirFilted);
hold on;
plot((ppgBeginning * 1 / fs - 1 / fs), ppglNirFilted(ppgBeginning),'rx','LineWidth',2);
title('PPG cycle beginnings');
xlabel('t [ms]');
ylabel('Blood pressure [mmHg]');

% RR-peaks
rrPeaks = pan_tompkins(ecgMvFiltered);

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

function out = toMillivolt(ECGsignal)
    out=3.3/8192*(ECGsignal-2048);
end

%convert the cuff pressure to mmHg
function out=toPressure(pressureSignal)
    out=(pressureSignal-175)/15.5;
end

 
function qrs=pan_tompkins(ppgf,fs)
%qrs=pan_tompkins(ekg)
%Csordas  Peter 07.03.26
%Pan Tompkins QRS detektor
% fs: optional - sampling freq in Hz - default is 1000
    if 1==nargin
        fs=1000;
    end
    fsh=fs/2;
    
    %Signal processing
    [b,a]=butter(5,5/fsh,'high'); ekgf=filtfilt(b,a,ppgf);                                   %Filtering 5-22 Hz
    [b,a]=butter(5,22/fsh); ekgf=filtfilt(b,a,ekgf);
    ekgd(1:length(ekgf),1)=0;                                                                 %Differentiate
    ekgd(3:end-2)= (-ekgf(1:end-4) -2*ekgf(2:end-3) +2*ekgf(4:end-1) + ekgf(5:end) );        
    %ekgp=abs(ekgd);                                                                         %Taking absolute value
    ekgp=ekgd.^2;
    N=round(150*fs/1000);
    ekgp=filter(ones(1,N)/N,1,ekgp);                                                      %150 ms Moving average

    
    %Peak detection
    tmp0=diff(ekgp);
    tmp1=ekgp.*([0; tmp0]>0).*([tmp0; 0]<0);
    peak_locations=find(tmp1); 

    %Ignore all peaks that precede or follow larger peaks by less than 200 ms.
    Tmin=round(200*fs/1000);
    tmp0=find(diff(peak_locations)<Tmin);
    while ~isempty(tmp0)
        tmp1=ekgp(peak_locations(tmp0))-ekgp(peak_locations(tmp0+1));
        peak_locations(tmp0+(tmp1>0))=0; peak_locations=nonzeros(peak_locations);
        tmp0=find(diff(peak_locations)<Tmin);
    end

%figure; plot(ekgp); linmark(peak_locations,'''k''');

%     tmp0=ekg(peak_locations);
%     idx=find((tmp0>500) & (tmp0<3500)); % A kiules nem csucs
%     peak_locations=peak_locations(idx);
    
    qrs=1;
    peak=ekgp(peak_locations);
    
    %Initial estimations
    SPKI=min(peak(1:10));
    NPKI=0;
    RRint(1:8)=360*fs/1000;
    
    for j=1:length(peak_locations)
        rr=mean(RRint);
        THRI=NPKI+0.35*(SPKI-NPKI); %Az eredeti algoritmusban 0.25 a kuszob
        rrnew=peak_locations(j)-qrs(end);
        %If the peak is larger than the detection threshold call it a QRS complex,or If no QRS has been detected within 1.5 R-to-R intervals, there was a peak that was larger than half the detection threshold, and the peak followed the preceding detection by at least 360 ms, classify that peak as a QRS complex.
        if peak(j)>THRI || (peak(j)>THRI/2 && rrnew>max(1.66*rr,360*fs/1000) )
                            
            %T szures, ha <360 az rr
            if rrnew<(360*fs/1000)
                tmp0=diff(ekgp(qrs(end):peak_locations(j)));
                if max(tmp0)<-0.8*min(tmp0)
                  NPKI=(7*NPKI+peak(j))/8;
                  continue;
                end
            end
            noise_cntr=0; %CsP
            qrs(end+1)=peak_locations(j);
            if (peak(j)>THRI)
                SPKI=(7*SPKI+peak(j))/8;
            else
                SPKI=(3*SPKI+peak(j))/4;
            end
            if rrnew>0.92*rr && rrnew<1.16*rr
                RRint=[RRint(2:end) rrnew];
            end                        
        else    %otherwise call it noise.
            NPKI=(7*NPKI+peak(j))/8;
            noise_cntr=noise_cntr+1; %CsP nagy QRS-nek vett csucs miatt SPKI elszallt? - fokozatos csokkentes
            if noise_cntr>3
                SPKI=7*SPKI/8;
            end
        end
    end
    qrs=qrs(3:end);
    
    %Csúcsba tolás
    for j=1:length(qrs)
        tmp0=max(1,qrs(j)-100); 
        tmp1=min(length(ekgf),qrs(j)+100);
        tmp1=ekgf(tmp0:tmp1); tmp1=find(tmp1==max(tmp1)); tmp1=tmp1(1);
        qrs(j)=tmp0+tmp1-1;
    end
    
    if size(ppgf,1)>1
        qrs=qrs'; %return the same type, as input
    end    
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