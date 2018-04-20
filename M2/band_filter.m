function filted=band_filter(type,order,raw)
% type: sz�r� t�pus, order: foksz�m, raw: nyers adat, filted: sz�rt adat
if strcmp(type,'butter')
[a,b]=butter(order,[0.09655172413 0.10049261083],'stop');
filted = filtfilt(a,b,raw);
end

if strcmp(type,'cheby')
[a,b]=cheby1(order,1,[0.09655172413 0.10049261083],'stop');
filted = filtfilt(a,b,raw);
end

end

