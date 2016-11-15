%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAKE TABLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 1: Preprocessing
load('/Users/PHUC_VO/Documents/ECE2200/FinalProject/sample.mat');
X = y(:,1);                     % use the first channel 
resampledSong = resample(X,8000,44100);     % Resample the signal at 8000Hz

% Step 2: Spectrogram 
Fs = 8000;                  % sampling rate of the signal
window = 64*10^(-3)*Fs;     % the length of chunks we take fft of
noverlap = 32*10^(-3)*Fs;   % number of samples overlap b/w adjacent chunks 
nfft = window;              % length od fft. In this case, nfft = window
[S,F,T] = spectrogram(resampledSong,window,noverlap,nfft,Fs);
log_S = log10(abs(S)+1);
%%% Check step 2:
%     figure 
%     imagesc(T,F,20*log10(abs(S)))
%     axis xy;
%     xlabel('Time (s)')
%     ylabel('Frequency (kHz)')
%     title('Spectrogram')
% 
%     colormap jet
%     c = colorbar;
%     set(c);
%     ylabel(c,'Power (dB)','FrontSize',14);

% Step 3: Feature Extraction (Spectrogram Local Peaks)
% log_S = [1 1 1 3;1 10 1 2;1 3 1 5;2 11 1 2];
localPeak = ones(size(log_S,1),size(log_S,2));  
for i = -4:4        % check magnitude with neighbor in the 9x9 grid (gs=9)
    for j = -4:4
        if (i~=0 || j~=0)
            CS = circshift(log_S,[j,i]);
            localPeak = and(localPeak,((log_S - CS)>0));
        end
    end
end
  

% Step 4: Thresholding
for i = 1:size(localPeak,1)         % if any peak < 0.5883 
    for j = 1:size(localPeak,2)     % then that peak will be removed
        if (localPeak(i,j)==1)
            if (log_S(i,j)<0.5883)
                localPeak(i,j)=0;
            end
        end
    end
end
nz = nnz(localPeak);

% Step 5: Constructing the table 
fanOut = 3;
deltaTL = 3;
deltaTU = 6;
deltaF = 9;



