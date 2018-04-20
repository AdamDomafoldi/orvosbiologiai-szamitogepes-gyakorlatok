function out = toMillivolt(ECGsignal)
    out = 3.3 / 8192 * (ECGsignal - 2048);
end