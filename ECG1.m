clc
clear all
close all

% ==============================================================================
%
%                               ECG Signal Processing
%
% ==============================================================================

% Load Noisy ECG data
load ECG_Lead2_900sec_Fs1000;

% ------------------------------------------------------------------------------
% Input Parameters.

noisySignal = Signal;
Fs = 1000;                               % Sampling Frequency.

% ------------------------------------------------------------------------------
% Low Pass Filtering Signals above 100Hz

% ------------------------------------------------------------------------------
% Filter has been designed using FDA tool with following parameters.

% Filter Type : Low pass FIR Equi ripple filter
% Pass band : 90 Hz
% Stop band : 110 Hz
% Pass band Ripple : 1 db
% Stop band Ripple : 80 db
% ------------------------------------------------------------------------------

% Load Filter Coefficients
load fc

lpFilteredSignal = filtfilt(lpfN,1,noisySignal);

% ------------------------------------------------------------------------------
% Displaying Low Pass filtered signal for 5 second duration

tStart = 50;                            % Starting time to display.
                                        % (Between 0 and 90 sec) 
n = (tStart*Fs)+1:(tStart + 5)*Fs;

subplot(2,1,1)
plot(n,noisySignal(n))
title('Noisy ECG Signal')
xlim([n(1),n(end)])

grid on
ylabel('Amplitude(mV)')
set(gca,'xTick',(tStart*Fs):Fs:(tStart + 5)*Fs,...
    'xTickLabel',tStart:1:tStart+5)


subplot(2,1,2)
plot(n,lpFilteredSignal(n))
title('Low Pass Filtered Signal')
xlim([n(1),n(end)])

xlabel('Time (sec) -->')
ylabel('Amplitude(mV)')

grid on
set(gca,'xTick',(tStart*Fs):Fs:(tStart + 5)*Fs,...
    'xTickLabel',tStart:1:tStart+5)

% Power Spectrum Estimation

nfft = 2^nextpow2(length(noisySignal));

Pnoisy = abs(fft(noisySignal,nfft)).^2/length(noisySignal)/Fs;
Hpsd = dspdata.psd(Pnoisy(1:length(Pnoisy)/2),'Fs',Fs);  

Plpf = abs(fft(lpFilteredSignal,nfft)).^2/length(lpFilteredSignal)/Fs;
Hpsd1 = dspdata.psd(Plpf(1:length(Plpf)/2),'Fs',Fs);

figure

subplot(2,1,1)
plot(Hpsd)
title('PSD of Noisy Signal')

subplot(2,1,2)
plot(Hpsd1)
title('PSD of Low Pass Filtered Signal')

% ------------------------------------------------------------------------------
% Notch Filtering to remove 60 Hz power line noise

% ------------------------------------------------------------------------------
% Filter has been designed using FDA tool with following parameters.

% Filter Type : Band Stop FIR Equiripple
% 
% Pass band Ripple : 1 db
% Stop band Ripple : 80 db

% ------------------------------------------------------------------------------

b =  [0.3576   -0.2385   -0.1831   -0.1484   -0.1297    0.8766   -0.1297 ...
     -0.1484   -0.1831  -0.2385    0.3576];

notchFilteredSignal = filtfilt(b,1,lpFilteredSignal);

% ------------------------------------------------------------------------------
% Displaying Notch filtered signal for 5 second duration

figure

subplot(2,1,1)
plot(n,noisySignal(n))
title('Noisy ECG Signal')
xlim([n(1),n(end)])

grid on
ylabel('Amplitude(mV)')
set(gca,'xTick',(tStart*Fs):Fs:(tStart + 5)*Fs,...
    'xTickLabel',tStart:1:tStart+5)


subplot(2,1,2)
plot(n,notchFilteredSignal(n))
title('Low Pass and Notch Filtered Signal')
xlim([n(1),n(end)])

xlabel('Time (sec) -->')
ylabel('Amplitude(mV)')

grid on
set(gca,'xTick',(tStart*Fs):Fs:(tStart + 5)*Fs,...
    'xTickLabel',tStart:1:tStart+5)

% Power Spectrum Estimation

Pnotch = abs(fft(notchFilteredSignal,nfft)).^2/length(notchFilteredSignal)/Fs;
Hpsd2 = dspdata.psd(Pnotch(1:length(Pnotch)/2),'Fs',Fs);

figure

subplot(2,1,1)
plot(Hpsd)
title('PSD of Noisy Signal')

subplot(2,1,2)
plot(Hpsd2)
title('PSD of Low Pass and Notch Filtered Signal')

% ------------------------------------------------------------------------------

% ==============================================================================
% Part 2 : Detection of R peaks
% ==============================================================================

t = (1:length(Signal))/Fs;
ECGsig = notchFilteredSignal;

%initialization of the matrices for storing R-peaks and location
R_peaks=zeros(1,length(t));
Index=zeros(1,length(t));

W = 400;    %window size represents number of samples in 0.4second (Fs = 1000)

for i = 1:W:length(t)
    
        [Rp,index] = max(ECGsig(i:i+W-1));  %finding the peak in each window.
        
        if Rp > 0.7     %storing the peak and location based on the criteria.
            
            R_peaks(i+index-1)=Rp;
            Index(i+index-1)=i+index-1;
            
        end
end

% Storing all non zero peaks and locations

nz = R_peaks ~= 0;

RPk = R_peaks(nz);
locs = Index(nz);

% chosing the higher peak if two peaks are detected in two adjacent windows
% with location difference shorter than W

for m = 1:length(locs)-1
    if locs(m+1) - locs(m) < W
        [Max_R, Index_R] = max(RPk(m:m+1));
        if Index_R == 1
            RPk(m+1) = 0;
            locs(m+1) = 0;
        else
            RPk(m) = 0;
            locs(m) = 0;
        end
    end
end

% removing all zero values from resultant matrix

nz = RPk ~= 0;

RPk = RPk(nz);
locs = locs(nz);

fprintf('Number of detected R Peaks = %d\n', length(RPk))
fprintf('Average Heart Rate over the duration = %.2f BPM\n',length(RPk)*60/900)

% ------------------------------------------------------------------------------
% Plotting Heart rate variability

t1 = t(locs);
ix = t1 > tStart & t1 <tStart + 20;
tn =  t1(ix);
hr = 60./(tn(2:end) - tn(1:end-1));
                       
n = (tStart*Fs)+1:(tStart + 20)*Fs;

figure

subplot(2,1,1)
plot(n,ECGsig(n))
hold on
plot(locs(ix),RPk(ix),'R*')
title('ECG Signal (Detected R Peaks)')
xlim([n(1),n(end)])
xlabel('Time(s)')

grid on
ylabel('Amplitude(mV)')
set(gca,'xTick',(tStart*Fs):2*Fs:(tStart + 20)*Fs,...
    'xTickLabel',tStart:2:tStart+20)

hold off

subplot(2,1,2)
plot([tStart tn tStart+20],[hr(1) hr(1) hr hr(end)],'-o','linewidth',2);
title('Heart rate variability vs Time')
grid on
xlabel('Time(s)')
ylabel('Heart rate')

% ------------------------------------------------------------------------------
% Area under QRS

qrsArea = zeros(1,length(RPk));

for i = 1:length(RPk)
    
    n0 = max(locs(i) - 50, 1) : locs(i) + 50;
    qrsArea(i) = sum(ECGsig(n0));
    
end

% Plotting

n = (tStart*Fs)+1:(tStart + 20)*Fs;

figure

subplot(2,1,1)
plot(n,ECGsig(n))
hold on
plot(locs(ix),RPk(ix),'R*')
iy = locs(ix);
for i = 1 : length(iy)
    plot(iy(i)-50:iy(i)+50,ECGsig(iy(i)-50:iy(i)+50),'m')
end
title('ECG Signal (QRS Area marked)')
xlim([n(1),n(end)])

xlabel('Time(s)')

grid on
ylabel('Amplitude(mV)')
set(gca,'xTick',(tStart*Fs):2*Fs:(tStart + 20)*Fs,...
    'xTickLabel',tStart:2:tStart+20)

q = qrsArea(ix);

subplot(2,1,2)
plot([tStart tn tStart+20],[q(1) q q(end)],'-o','linewidth',2);
title('QRS Area under curve Vs Time')
grid on
xlabel('Time(s)')
ylabel('QRS Area')

% ==============================================================================

    
    
    










