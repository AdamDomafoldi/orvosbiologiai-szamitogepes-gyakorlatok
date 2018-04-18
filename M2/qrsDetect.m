function [ tR ] = qrsDetect( rawECG )
%QRSDETECT Summary of this function goes here
%   Detailed explanation goes here
close all;
fs=1000;
t=(0:length(rawECG)-1)/fs;
fl=6/500;
fh=16/500;
[fBnum fBden]=butter(3,[fl fh],'bandpass');
figure();
freqz(fBnum, fBden); %frekvencia valasz abrazolasa
ecgbp=filter(fBnum,fBden,rawECG);
figure();
plot(t,ecgbp);
hold on;
derivate=zeros(1,length(ecgbp));
derivate(1:length(ecgbp)-1)=diff(ecgbp);
plot(t,derivate,'k');
d2=derivate.*derivate;    
plot(t,d2,'c');
sampleInt=round(0.2*fs); %mert 200 ms a qrs hossza és ott vannak nagy meredeksegek;
sInt=zeros(1,length(d2));
for i=sampleInt:length(d2)-sampleInt;
   sInt(i)= sum(d2(i-round(sampleInt/2):i+round(sampleInt/2))); 
end
plot(t,sInt,'r');

%% megkereses az integralt jelen a csucsokat
tre=0.005;
w=0.15*fs;
j=1;
for i=1:length(sInt)
   if( sInt(i)>tre && sInt(i)==max(sInt(i-w:i+w)) )
   tR(j)=t(i);
   j=j+1;
   end
       
end
plot(tR,0.2,'rx');
legend('filt ecg','derivate','derivate square','integralt jel','talalt QRS');


end

