close all; clear all;
aa=hhmbinread('brigitta_relaxed.hhm');
ppg=aa.ppgl_nir;
ecg=toMillivolt(aa.ecg1);

%Tényleges mintavételi fekvencia
[f , ekg_spect] = toSpect(ecg);
[M,I] = max(ekg_spect);
real_50=f(I);
fs=(50/real_50)*1000;

ecg=ecg(5*fs:end); %kezdo tranziens levagasa
ppg=ppg(5*fs:end);
t=(1:length(ecg))/fs;
figure();
plot(t,ecg);


%% hrv
tR=pan_tompkins(ecg,1000); %ms-ben van most
tRR=diff(tR);
MSD=mean(tRR);
SDNN=std(tRR);
RMSSD=sqrt(mean(tRR.*tRR));
j=1;
for i=2:length(tRR)
   if(tRR(i)-tRR(i-1)>50)
      j=j+1; 
   end
end
pNN50=j/length(tRR);
poinx=tRR(1:end-1);
poiny=tRR(2:end);
figure();
temp=linspace(min(tRR),max(tRR),1000);
plot(poinx,poiny,'b.',temp,temp,'m', temp,temp+370,temp,temp-281);
title('poincare diagramm');
xlabel('tRR(i) [ms]'); ylabel('tRR(i+1) [ms])');

%% 2. példa
close all;
figure();
%50/500=0.1 az 50Hz
[f50num f50den]=butter(3,[0.09 0.11],'stop');

ecg50=filtfilt(f50num,f50den,ecg);
[fn fd]=butter(3,[0.5/500 35/500],'bandpass'); 
ecg50=filtfilt(fn,fd,ecg50);
plot(t,ecg50,'b',tR/fs,ecg50(tR),'rx'); hold on;
% T hullámok keresése
% R-R közötti lokális maximum az a T
w=0.1*fs;
tR(1)=[];
for i=1:length(tR)-1
   [maxtmp, maxIndex]=max(ecg50( tR(i)+w : tR(i)+tRR(i)/2) );
   tT(i)=tR(i)+w+maxIndex;
   [mintmp, minIndex]=min(ecg50(tT(i):tT(i)+150));
   tT(i)=tT(i)+minIndex; % a T hullám végét keressuk meg így
%    fele(i)=tR(i+1)-40;%tR(i+1)-round( tRR(i)/2 );
   [maxtmp, maxIndex]=max(ecg50( tR(i)+tRR(i)/2 : tR(i+1)-40 ));
   tP(i)=tR(i)+round(tRR(i)/2)+maxIndex;
   [mintmp, minIndex]=min(ecg50( tP(i)-60 : tP(i) ) );
   minIndex=60-minIndex;
   tP(i)=tP(i)-minIndex;
   
   %q hullám megkeresese
   [mintmp minIndex]=min(ecg50( tR(i)-60 : tR(i) ) );
   minIndex=60-minIndex;
   tQ(i)=tR(i)-minIndex;
   
   
end
if(tP(1)>tQ(1))
    tQ(1)=[];
end
if( tP(1)>tR(1) )
    tR(1)=[];
end
if( tP(1)>tT(1) )
    tT(1)=[];
end
tP=tP(1:end-2);
tQ=tQ(1:length(tP));
tR=tR(1:length(tP));
tT=tT(1:length(tP));


close all;
plot(t,ecg50,'b'); %always
hold on;
plot(tR/fs,ecg50(tR),'rx');
plot(tT/fs,ecg50(tT),'g*');
% plot(fele/fs,ecg50(fele),'k*');
plot(tP/fs,ecg50(tP),'rd');
plot(tQ/fs,ecg50(tQ),'ms');
title('A talalt pontok');
xlabel('t [s]'); ylabel('ECG [mV]'); legend('EKG');
legend('ecg','R','T','P','Q');

PQ=(tQ-tP);
QT=tT-tQ;
RR=diff(tR);
meanPQ=mean(PQ); meanQT=mean(QT); meanRR=mean(RR); 
stdPQ=std(PQ); stdQT=std(QT); stdRR=std(RR); 
%% 3 4 példa
close all;

ppg=ppg(8500:14000);
ecg50=ecg50(8500:14000);
clear tR;
tR=pan_tompkins(ecg50);
tR=tR';


[blalbalab tPPG]=findpeaks(-ppg,'MinPeakDistance',600,'MinPeakHeight',-2000);

t=(0:length(ppg)-1)/fs;
subplot(2,1,1);
plot(t,ppg,'r',tPPG/fs,ppg(tPPG),'bx');
title('PPG jel');
xlabel('t [s]'); ylabel('Amplitude'); 
hold on;
subplot(2,1,2);
plot(t,ecg50,'b',tR/fs,ecg50(tR), 'rx');
title('ECG');
xlabel('t [s]'); ylabel('ECG [mV]');

% %szivciklus atlagos hossza
% hossz=mean(diff(tPPG));
% kul=mean(tPPG-tR);
% Speed=0.93/(mean(tPPG-tR)/fs);

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












