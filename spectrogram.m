function [FData,freq,time] = spectrogram(data,Fs,timeWindow,timeOverlap,maxf,db)
% This function was made by Mike, NOT a Matlab built in function
%
% [FData,freq,time] = spectrogram(data,Fs,timeWindow,timeOverlap,maxf,db)
% calling spectrogram with no output produces plot of 
% signal's frequency contents over time
% Inputs data, Fs, timeWindow, timeOverlap, and maxf must be specified
% if db is not set the plot will be in absolute scale, set db to 1 
% for log scale plot
%
% outputs FData is a matrix where column are the frequency content amplitudes
% at fft times. freq and time are the frequency and time vectors respectively
% 
% Try not to make timeWindow ludicriously large, try to limit it to a fourth
% the length of the data vector

if timeOverlap >= timeWindow
    error('overlap must be less than time window')
end

if maxf > Fs/2, error('maxf must be less than half of sampling frequency'); end

if size(data,1) == 1 ,tempdata = data'; else tempdata = data; end

if ~exist('db','var')
    db = 0;
else
    if (db ~= 1)&& (db ~= 0)
        error('1 for db scale, 0 for absolute scale')
    end
end
    
%moving fft
sampWindow = floor(timeWindow*Fs);
sampOverlap = floor(timeOverlap*Fs);
increment = sampWindow - sampOverlap;

FData = [];
time = zeros(1,1);

index = 1;
while index + sampWindow <= length(data)
    column = 2*abs(fft(tempdata(index:index+sampWindow-1)))/sampWindow;
    if db == 0
        FData = [FData, column];
    else
        FData = [FData, 20*log10(column)];
    end
    index = index + increment;
    time = [time, time(end)+increment/Fs];
end

%truncate matrix to view relevant bandwidth
freqResolution = Fs/sampWindow;
FData = FData(1:floor(maxf/freqResolution),:);

%time and frequency vectors
freq = (0:Fs/sampWindow:Fs);
freq = freq(1:ceil(maxf/freqResolution));
time = time(1:end-1);

%plot
%mesh(time,freq,FData) 
surf(time,freq,FData,'EdgeColor','none','FaceColor','interp')
xlabel('time [s]')
ylabel('frequency [Hz]')
view(119.3,42.8)

%set colormap
cMap = parula(256);
if db == 0                      %colormap for absolute scale
    dataMax = 0.25;
    dataMin = 0;
    centerPoint = 0.07;
    scalingIntensity = 4;
    zlabel('amplitude [g]')
else                            %colormap for db scale
    dataMax = -10;
    dataMin = -100;
    centerPoint = -50;
    scalingIntensity = 3;
    zlabel('amplitude [db]')
end
x = 1:length(cMap);
x = x - (centerPoint-dataMin)*length(x)/(dataMax-dataMin);
x = scalingIntensity * x/max(abs(x));
x = sign(x).* exp(abs(x));
x = x - min(x); x = x*511/max(x)+1;
newMap = interp1(x, cMap, 1:512);
colormap(newMap)

end
