%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAKE TABLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc;
% Step 1: Preprocessing
load('/Users/PHUC_VO/Documents/ECE2200/FinalProject/sample.mat');
X = y(:,1);                 % use the first channel 
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
gs = 9;
for i = -(gs-1)/2:(gs-1)/2          % check magnitude with neighbor 
    for j = -(gs-1)/2:(gs-1)/2      % in the 9x9 grid (gs=9)
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
[I,J] = ind2sub(size(localPeak),find(localPeak));% coordinate of peaks in localMatrix
peakIn = [I'; J'];                               % storage in matrix peakIn

count = 0 ;
F1 = 0; F2 = 0; T1=0; T2=0;

for peak = 1:size(peakIn,2)
    F_start = peakIn(1,peak)-deltaF;            % Create 4 points of 
    F_end = peakIn(1,peak)+deltaF;              % sub matrix
    T_start = peakIn(2,peak)+deltaTL;           % (region to get 3 peaks
    T_end = peakIn(2,peak)+deltaTU;             % to compare)
    
    if (F_start <= 0)       % Re-assign value for bottom-left point of
        F_start = 1;        % sub matrix if starting point <= 0
    end
    if (F_end > 257)        % Re-assign value for top-left point of 
        F_end = 257;        % sub matrix if ending point > 257
    end 
    if(T_end<=936)          % Check the right-side of submatrix is
        count = count+1;    % not exceed the localMatrix dimension 
        
        subMatrix = localPeak(F_start:F_end,T_start:T_end); % create submatrix
        
        % find and save peaks in submatrix 
        % storage peaks "coordinates in submatrix" into to matrix A
        [row,col] = ind2sub(size(subMatrix),find(subMatrix,fanOut));
        A = [row';col'];
    
        if(nnz(A)~=0)   % if there are peaks in submatrix
            
         % convert "coordinate in submatrix" to get coordinates in localPeak
         if (peakIn(1,peak)>deltaF) % only convert row when bottom-left point
           row = row+F_start-1;     % of submatrix < deltaF (9)
         end 
         col = col+T_start-1;       % convert column coordinate
         if((row>257) | peakIn(1,peak)>257 | (col>936) | peakIn(2,peak)>936)
             count = count+1;
         end
            % If there is one peak found in submatrix
            F1 = [F1 F(peakIn(1,peak))];
            F2 = [F2;F(row)];
            T1 = [T1 T(peakIn(2,peak))];
            T2 = [T2 T(col)];
            
         % If there are 2 peaks found in submatrix
         if(nnz(A)/2==2)
            F1 = [F1 F(peakIn(1,peak))];
            T1 = [T1 T(peakIn(2,peak))];
         end 
         
         % If there are 3 peaks found in submatrix
         else if(nnz(A)/2==3)
            F1 = [F1 F(peakIn(1,peak))];
            F1 = [F1 F(peakIn(1,peak))];
            T1 = [T1 T(peakIn(2,peak))];
            T1 = [T1 T(peakIn(2,peak))];
         end 
        end
    end
end

% Step 6: Final Function
table = [F1;F2';T1;T2-T1];
table = table';
table = table(2:size(table,1),:);
count 




