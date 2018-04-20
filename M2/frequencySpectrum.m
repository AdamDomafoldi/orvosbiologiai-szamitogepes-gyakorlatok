function frequencySpectrum(ecg,titleOfDiagram)
    % @ecg: hhm object's property
    % @titleOfDiagram: title of the diagram
    fs=1.0153e+03; % sampled at this frequency [Hz]   
    Y = fft(ecg); % fast fourier transformation
    L=length(ecg); % signal length
    signal = abs(Y/L); % absolute value of frequency
    spektrum = signal(1:L/2+1); % cut the signal into half
    f = fs *(0:(L/2))/L; % frequency axsis [Hz] 
    % visualize half signal
    figure();
    plot(f,spektrum) 
    title(titleOfDiagram)
    xlabel('f [Hz]')
    ylabel('amplification')
    xlim([0 500]) % X axis
end