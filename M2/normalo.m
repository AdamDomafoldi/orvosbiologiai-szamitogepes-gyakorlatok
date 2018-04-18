function [ avgECG ] = normalo( rawECG )
%NORMALO Summary of this function goes here
%   Detailed explanation goes here
fs=1000;
t=(0:length(rawECG)-1)/fs;
tR=qrsDetect(rawECG);
close all;
firstPeriod=rawECG(tR(1)*fs:tR(2)*fs);
t1=(0:length(firstPeriod)-1)/fs;
figure;
plot(t1,firstPeriod);
N=length(firstPeriod);
avgECG=firstPeriod;
for i=2:length(tR)-1
    actPeriod=rawECG(tR(i)*fs:tR(i+1));
    [num,den]=rat(N/length(actPeriod));
    
    avgECG=(avgECG+normalized)/2;
end

end

