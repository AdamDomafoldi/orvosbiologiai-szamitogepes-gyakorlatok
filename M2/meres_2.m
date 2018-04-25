clear all;

%% 1.feladat

struct=hhmbinread('brigitta_relaxed.hhm'); % adatok beolvasása
ekg=struct.ecg1;
ekg=ekg(10000:size(ekg,1));
ekg=toMillivolt(ekg);


%Tényleges mintavételi fekvencia
[f , ekg_spect] = toSpect(ekg);
[M,I] = max(ekg_spect);
real_50=f(I);
fs=(50/real_50)*1000;


ekg_time=((size(ekg,1))-1)/fs;
t=0:1/fs:ekg_time;
plot(t,ekg);

[B1,A1]=butter(2,0.0005,'low'); %2.fokú aluláteresztõ Butterworth szûrõ
ekgalap = filtfilt(B1,A1,ekg);%szûrt jel

figure
hold on
plot(t,ekg,'b');
plot(t,ekgalap,'r','Linewidth',3);
hold off

min_ekg = min(ekgalap);
max_ekg = max(ekgalap);
alap_inter = max_ekg-min_ekg;
disp('Alapjel vándorlás intervalluma:');
disp(alap_inter);


%% 2. feladat

% EKG jel frekvencia tartományban

frequencySpectrum(ekg,'Nyers jel frekvencia tartomány');

%Szûrés
ekg_50_szurt=band_filter('butter',2,ekg);
frequencySpectrum(ekg_50_szurt,'Szûrt frekvencia tartomány');

% Zaj
[B2,A2]=butter(2,[0.09655172413 0.10049261083],'bandpass'); %50Hz-es zaj
zaj=filtfilt(B2,A2,ekg);
frequencySpectrum(zaj,' Zaj frekvencia tartomány');

% Jel-zaj viszony
jel = ekg-zaj;
jel_szurt = ekg_50_szurt-zaj;
jel_zaj = snr(jel,zaj); %jel-zaj arány szûrés elõtt

% Vizsgálat: szûrés után maradt-e zaj
zaj_szurt=band_filter('butter',2,zaj);
frequencySpectrum(zaj_szurt,' Szûrt zaj frekvencia tartomány'); %Maradt zaj, de nagyon kis mértékben

jel_zaj_szurt = snr(jel_szurt,zaj); % Szûrés utáni jel-zaj arány


%% 3.feladat

[B3,A3]=butter(2,0.06896551724,'low'); %2.fokú aluláteresztõ Butterworth szûrõ
ekg_35 = filtfilt(B3,A3,ekg);%szûrt jel

figure

hold on
plot(t,ekg,'b');
plot(t,ekg_35,'r');
hold off

%% 4.feladat

[B4,A4]=butter(2,[0.01 0.03],'bandpass');
ekg_QRS_szurt=filtfilt(B4,A4,ekg);

[R_ertekek , R_hely]=findpeaks(ekg_QRS_szurt,'MinPeakDistance',500);

figure 
hold on
plot(t,ekg_QRS_szurt,'b');
scatter(R_hely/fs,R_ertekek,'r');
hold off

%% 5. feladat

% ez nem volt jo :/ rip algoritmus
% dist = (148);
% 
%  for i=1:length(R_hely)
%     if i+2 == length(R_hely)
%         break 
%     end
%     dist(i)=((R_hely(i+2)-R_hely(i+1))/2+R_hely(i+1)) - ((R_hely(i+1)-R_hely(i))/2+R_hely(i));
%  end
 
% 5. Feladat:
period=500; % erre az alapperiódusra normálunk, 500 minta
RR_dist = R_hely(2:end) - R_hely(1:end-1); % különbözõ RR távok
RR_interpol = zeros(length(RR_dist), period);
for i=1:length(RR_dist)
step = (RR_dist(i)-1) / (period-1); % lépésköz
val = interp1([1:RR_dist(i)+1], ekg(R_hely(i):R_hely(i+1)),[1:step:RR_dist(i)]);
RR_interpol(i,:)=val;
end
atlagolt = mean(RR_interpol);
zajos = RR_interpol(1,:);
figure;
plot([zajos(251:500) zajos(1:250)]);
% sorrend csere, hogy ne csúcstól csúcsig nézzük
title('Zajos alapperiódus');

figure;
plot([atlagolt(251:500) atlagolt(1:250)]);
title('Átlagolt EKG');

% EKG jel frekvencia tartományban

frequencySpectrum(ekg,'Nyers jel frekvencia tartomány');

%Szûrés
frequencySpectrum(atlagolt,'Szûrt frekvencia tartomány');

% Zaj
[B2,A2]=butter(2,[0.09655172413 0.10049261083],'bandpass'); %50Hz-es zaj
zaj=filtfilt(B2,A2,atlagolt);
frequencySpectrum(zaj,' Zaj frekvencia tartomány');

% Jel-zaj viszony
jel_szurt_uj = atlagolt-zaj;
jel_zaj_szurt_uj = snr(jel_szurt_uj,zaj); % Szûrés utáni jel-zaj arány
