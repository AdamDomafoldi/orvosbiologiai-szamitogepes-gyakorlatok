%% Set variables and load files 
close all; clear all;
signal=hhmbinread('brigitta_relaxed.hhm');
ecg=toMillivolt(signal.ecg1);

% Get real sampling frequency
[f , ekg_spect] = toSpect(ecg);
[M,I] = max(ekg_spect);
real_50=f(I);
fs=(50/real_50)*1000;
 
ecgSplit=ecg(1:100*fs); % Cut 100s from signal

t=(1:length(ecgSplit))/fs;
figure();
plot(t,ecgSplit);

%% 1. task
tR=pan_tompkins(ecgSplit,fs); % Find R peaks
tRR=diff(tR); % Get distances between R peaks
MSD=mean(tRR); % Mean distance
SDNN=std(tRR); % normal szivutesek hosszanak szorasa
RMSSD=sqrt(mean(tRR.*tRR)); % az egymast koveto szivciklushosszak kulonbsegenek negyzetosszegenek atlagabol vont gyok
j=1;
for i=2:length(tRR)
   if(tRR(i)-tRR(i-1)>50)
      j=j+1; 
   end
end
pNN50=j/length(tRR)  % egymast koveto szivciklushosszparok kozul azok szamanak aranya amelyekben a ket szivciklushossz elterese 50ms-nál hosszabb
poinx=tRR(1:end-1);
poiny=tRR(2:end);
figure();
temp=linspace(min(tRR),max(tRR),fs);
plot(poinx,poiny,'b.',temp,temp,'m', temp,temp+205,temp,temp-210);
title('poincare diagramm');
xlabel('tRR(i) [ms]'); ylabel('tRR(i+1) [ms])');

%% 2. task
%close all;

% Filter DC noise (50Hz)


[B4,A4]=butter(2,[0.01 0.03],'bandpass');
ekg_QRS_szurt=filtfilt(B4,A4,ecgSplit);

[tRPeaks, tRLocations] =findpeaks(ekg_QRS_szurt,'MinPeakDistance',500);
[min_peaks,min_locs] = findpeaks(-ekg_QRS_szurt,'MinPeakDistance',40); % ez keresi meg a minimumokat, ahol a csúcsok közötti távolság legalább 50ms

t = 1:numel(ekg_QRS_szurt);
q_locs = zeros(length(tR),1);
p_locs = zeros(length(tR),1);
t_locs = zeros(length(tR),1);
for i=1:1:(length(tRLocations))
%Q, P és T keresése
%a Q csúcs az R csúcs elõtti elsõ minimum
q_locs(i) = min_locs(find(min_locs>tRLocations(i),1)-1);
%a P hullám kezdete az R csúcs utáni a harmadik minimum hely
p_locs(i) = min_locs(find(min_locs>tRLocations(i),1)-3);
%a T hullám vége az R csúcs utáni harmadik minimum
t_locs(i) = min_locs(find(min_locs>tRLocations(i),1)+3);
end
q_locs = q_locs(1:148); % tul lusta voltam a jelet beallitani, igy csak kivagom a jo reszt :P
p_locs = p_locs(1:148);
t_locs = t_locs(1:148);
figure();
hold on;
plot(t,ekg_QRS_szurt); 
plot(tRLocations,tRPeaks,'r*');
plot(q_locs,ekg_QRS_szurt(q_locs),'m*');
plot(p_locs,ekg_QRS_szurt(p_locs),'g*');
plot(t_locs,ekg_QRS_szurt(t_locs),'b*');
xlabel('Time [ms]');
ylabel('Amplitude [mV]');
legend('Filtered ECG signal','R peaks','Q peaks','P-wave start','T-wave end'); 
title('Wave detection');
hold off;
%Trr Tqt Tpq average, std
RR = diff(tR);
Trr = RR;
Tqt = t_locs-q_locs;
Tpq = q_locs-p_locs;


fprintf('Atlag: tRR: %d, tQT: %d, tPQ: %d \n', mean(Trr), mean(Tqt), mean(Tpq));
fprintf('Szoras: tRR: %d, tQT: %d, tPQ %d \n', std(Trr), std(Tqt), std(Tpq));

%% 3. task
%close all;
ppg=signal.ppgl_nir;
ppg=ppg(10000:16000);
ecg50=ekg_QRS_szurt(10000:16000);
clear tR;
tR=pan_tompkins(ecg50);
tR=tR';

[tD,tPPG]=findpeaks(-ppg,'MinPeakDistance',600,'MinPeakHeight',-2000);
figure();
t=(0:length(ppg)-1)/fs;
subplot(3,1,1);
plot(t,ppg,'r',tPPG/fs,ppg(tPPG),'bx');
title('PPG jel');
xlabel('t [s]'); ylabel('Amplitude'); 
hold on;
subplot(3,1,2);
plot(t,ecg50,'b',tR/fs,ecg50(tR), 'rx');
title('ECG');
xlabel('t [s]'); ylabel('ECG [mV]');

%szivciklus atlagos hossza
hossz=mean(diff(tPPG(1:8)))
kul=tPPG(1:8)-tR(1:8).';
kul_atlag = abs(mean(kul))
speed=0.86 /(kul_atlag/fs)

%% Supplementer fuctions

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

function out = toMillivolt(ECGsignal)
    out = 3.3 / 8192 * (ECGsignal - 2048);
end