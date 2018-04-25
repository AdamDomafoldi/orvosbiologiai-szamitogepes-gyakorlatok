function [ cycles ] = CycleDetect( be )
close all
fs=1000;
be=be(25*fs:end);
t=(0:length(be)-1/fs);



hold on
plot(t,be);


sampleInt=round(0.8*fs);
sInt=zeros(1,length(be));

plot(t,sInt,'r');

for i=sampleInt:length(be)-sampleInt;
   sInt(i)= sum(be(i-round(sampleInt/2):i+round(sampleInt/2))); 
end
    
max(sInt)/0.5;
w=0.5*fs;
j=1;

for i=fs+100:length(sInt)-2*fs
   tre=max(sInt(i-fs/2:i+fs/2))/2; 
   if( sInt(i)>tre  && be(i)==min(be(i-w:+i+w)) );
   tR(j)=t(i);
   j=j+1;
   end
   
end

plot(tR,0.8,'rx');
xlabel('time [s]'); ylim([-500 4000]);

end
