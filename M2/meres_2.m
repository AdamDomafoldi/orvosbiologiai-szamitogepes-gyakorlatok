clear all;

%% 1.feladat

struct=hhmbinread('brigitta_relaxed.hhm'); % adatok beolvas�sa
ekg=struct.ecg1;
ekg=ekg(10000:size(ekg,1));
ekg=toMillivolt(ekg);


%T�nyleges mintav�teli fekvencia
[f , ekg_spect] = toSpect(ekg);
[M,I] = max(ekg_spect);
real_50=f(I);
fs=(50/real_50)*1000;


ekg_time=((size(ekg,1))-1)/fs;
t=0:1/fs:ekg_time;
plot(t,ekg);

[B1,A1]=butter(2,0.0005,'low'); %2.fok� alul�tereszt� Butterworth sz�r�
ekgalap = filtfilt(B1,A1,ekg);%sz�rt jel

figure
hold on
plot(t,ekg,'b');
plot(t,ekgalap,'r','Linewidth',3);
hold off

min_ekg = min(ekgalap);
max_ekg = max(ekgalap);
alap_inter = max_ekg-min_ekg;
disp('Alapjel v�ndorl�s intervalluma:');
disp(alap_inter);


%% 2. feladat

% EKG jel frekvencia tartom�nyban

frequencySpectrum(ekg,'Nyers jel frekvencia tartom�ny');

%Sz�r�s
ekg_50_szurt=band_filter('butter',2,ekg);
frequencySpectrum(ekg_50_szurt,'Sz�rt frekvencia tartom�ny');

% Zaj
[B2,A2]=butter(2,[0.09655172413 0.10049261083],'bandpass'); %50Hz-es zaj
zaj=filtfilt(B2,A2,ekg);
frequencySpectrum(zaj,' Zaj frekvencia tartom�ny');

% Jel-zaj viszony
jel = ekg-zaj;
jel_szurt = ekg_50_szurt-zaj;
jel_zaj = snr(jel,zaj); %jel-zaj ar�ny sz�r�s el�tt

% Vizsg�lat: sz�r�s ut�n maradt-e zaj
zaj_szurt=band_filter('butter',2,zaj);
frequencySpectrum(zaj_szurt,' Sz�rt zaj frekvencia tartom�ny'); %Maradt zaj, de nagyon kis m�rt�kben

jel_zaj_szurt = snr(jel_szurt,zaj); % Sz�r�s ut�ni jel-zaj ar�ny


%% 3.feladat

[B3,A3]=butter(2,0.06896551724,'low'); %2.fok� alul�tereszt� Butterworth sz�r�
ekg_35 = filtfilt(B3,A3,ekg);%sz�rt jel

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
period=500; % erre az alapperi�dusra norm�lunk, 500 minta
RR_dist = R_hely(2:end) - R_hely(1:end-1); % k�l�nb�z� RR t�vok
RR_interpol = zeros(length(RR_dist), period);
for i=1:length(RR_dist)
step = (RR_dist(i)-1) / (period-1); % l�p�sk�z
val = interp1([1:RR_dist(i)+1], ekg(R_hely(i):R_hely(i+1)),[1:step:RR_dist(i)]);
RR_interpol(i,:)=val;
end
atlagolt = mean(RR_interpol);
zajos = RR_interpol(1,:);
figure;
plot([zajos(251:500) zajos(1:250)]);
% sorrend csere, hogy ne cs�cst�l cs�csig n�zz�k
title('Zajos alapperi�dus');

figure;
plot([atlagolt(251:500) atlagolt(1:250)]);
title('�tlagolt EKG');

% EKG jel frekvencia tartom�nyban

frequencySpectrum(ekg,'Nyers jel frekvencia tartom�ny');

%Sz�r�s
frequencySpectrum(atlagolt,'Sz�rt frekvencia tartom�ny');

% Zaj
[B2,A2]=butter(2,[0.09655172413 0.10049261083],'bandpass'); %50Hz-es zaj
zaj=filtfilt(B2,A2,atlagolt);
frequencySpectrum(zaj,' Zaj frekvencia tartom�ny');

% Jel-zaj viszony
jel_szurt_uj = atlagolt-zaj;
jel_zaj_szurt_uj = snr(jel_szurt_uj,zaj); % Sz�r�s ut�ni jel-zaj ar�ny
