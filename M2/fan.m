function [ fanX fanY ] = fan( ecg,e )
%FAN Summary of this function goes here
%   Detailed explanation goes here
j=2;
% e=0.2;
i=1;
t=(0:length(ecg)-1)/1000;
savedX(1)=t(1);
savedY(1)=ecg(1);
fanY(1)=ecg(1);
fanX(1)=t(1);

slopeU=ecg(i+1)+e-savedX;
slopeL=ecg(i+1)-e-savedX;
    
while i <= length(ecg)-10
    k=2;
    limitH=slopeU*k;
    limitL=slopeL*k;
    while( ecg(i+k)<savedY+slopeU*k && ecg(i+k)>savedY+slopeL*k && i+k<(length(ecg)-10)) %benn van az intervallumban eldobjuk
        slopeUnew=(ecg(i+k)+e-savedX)/k;
        slopeLnew=(ecg(i+k)-e-savedX)/k;
        
        slopeU=min([slopeU slopeUnew]);
        slopeL=max([slopeL slopeLnew]); %uj meredeksegek letrehozasa
    limitH=slopeU*k;
    limitL=slopeL*k;
    
        k=k+1;
    end
    i=i+k;
    savedY=ecg(i);
    slopeU=ecg(i+1)+e-savedX;
    slopeL=ecg(i+1)-e-savedX; %pontot mentettunk uj meredeksegre van szuksegunk

    
    fanY(j)=savedY;
    fanX(j)=t(i);
    
    j=j+1;
end
close all;
figure();
plot(t,ecg,'b',fanX,fanY,'rx');
legend('eredeti','fan tomoritett');


end

