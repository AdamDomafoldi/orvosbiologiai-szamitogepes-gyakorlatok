%HHM file reader
% Input: name of HHM file (with full path and extension)
% Output: structure with fields:
%  ecg1: ECG lead Einthoven I.
%  ecg2: ECG lead Einthoven II.
%  press: cuff pressure captured from analog/digital converter. Offset is
%     not compensated. 15.5 LSB/mmHg
%  ppgl_red: PPG left red
%  ppgl_red_dc: PPG left red without DC level filtering (smaller gain)
%  ppgr_nir: PPG right near infra red
%  ppgl_nir_dc: PPG left nir without DC level filtering (smaller gain)
%  ppgl_nir: PPG left near infra red
%
% Every channel is sampled with 1 kHz
%Csordás Péter 2006
function o=hhmbinread(filename)
fid=fopen(filename,'r','b');
if fid==-1 error('File not found'); end

buffer=fread(fid,'bit12=>double'); 
idx=find(buffer<0);
if(~isempty(idx))
    buffer(idx)=buffer(idx)+4096; %overflow correction
end

o.ecg1=buffer(1:8:end);
o.ecg2=buffer(2:8:end);
o.press=buffer(3:8:end);
o.ppgl_red=buffer(4:8:end);
o.ppgl_red_dc=buffer(5:8:end);
o.ppgr_nir=buffer(6:8:end);
o.ppgl_nir_dc=buffer(7:8:end);
o.ppgl_nir=buffer(8:8:end);


