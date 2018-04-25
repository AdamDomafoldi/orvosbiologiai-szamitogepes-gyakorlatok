function [f,spektrum] = toSpect(ecg)
    % @ecg: hhm object's property
    % @titleOfDiagram: title of the diagram
    fs=1000; % sampled at this frequency [Hz]   
    Y = fft(ecg); % fast fourier transformation
    L=length(ecg); % signal length
    signal = abs(Y/L); % absolute value of frequency
    spektrum = signal(1:L/2+1); % cut the signal into half
    f = fs *(0:(L/2))/L; % frequency axsis [Hz] 
end