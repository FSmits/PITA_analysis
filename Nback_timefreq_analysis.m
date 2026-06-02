
clear
close all

% open EEGlab
cd('/Users/fsmits2/Downloads/eeglab2021.0')
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;


% Define EEG variables EEG.srate and EEG.time (same for every subject after preproc)
srate = 256; % in Hz

% Time vector that is mean-centered (-1 to +1)
time = (0:2 * srate)/srate;
time = time - mean(time);

% frequency parameters
min_freq =  2; % in Hz
max_freq = 30; % in Hz
num_freq = 56; % in count
% frequencies vector log spaced (because you want to focus on lower frequencies)
frex = logspace(log10(min_freq),log10(max_freq),num_freq);

% set a few different wavelet widths ("number of cycles" parameter)
% (also logarithmic spaced, like the frex)
num_cycles = logspace(log10(3), log10(15), num_freq); %doing the convolution with different wavelets per frequency, so lower number of cycles for lower freqs, and higher no. cycles for higher freqs.
% %Also possible to specify in time dimension (als fmhm):
% fwhms = linspace(.5,.2,numfrex);

% For baseline correction (decibel):
% specify baseline periods for dB-normalization
baseline_windows = [ -350 -100];

% Subject nrs in vector
Subj = ones(1,69)*110000; 
Subji = [1:2 4:6 8 10:12 14:27 29:31 34:50 52:55 57:67 69:79]; % dropouts=110007, 110013, 110032, 110051, 110068; misperformers=110028, 110033, 110056; no post-assess EEG: 110003, 110009
Subj = Subj + Subji;

for Subjecti = 1:length(Subj)
    
    CurrSubj = Subj(Subjecti);
    
    % Load data:
    LoadName = [num2str(CurrSubj) 'N-Back' 'Nameting' '-preproc-rmeyeblinkICA' '.set']
    EEG = pop_loadset( 'filename',LoadName,...
        'filepath','/Volumes/HER-1/onderzoeksarchief/19-452_CONTROL_BS/F_DataAnalysis/NBack EEG/Raw data sets NBack');

    %----------------------------------------------------------------------
    %            Time-freq analysis based on course Mike X Cohen
    %----------------------------------------------------------------------
    
    % Select the right epochs by finding which trial belongs to which condition:
    Epochs      = find( str2double({EEG.event.type}) > 3000 );
      
    % N's of convolution
    nData = prod(size(EEG.data(1,:,:)));
    nKern = length(time);
    nConv = nData + nKern - 1;
    halfK = floor(nKern/2);
    
    %%% Do the fft for data in all channels at once in matrix:
    % Use data in double precision for better computation precision
    data = double( EEG.data );
    % First reshape the data to be supertrials - all trials concatenated:
    dataR = reshape(data, EEG.nbchan, []);
    % FFT on data (2nd dimension - timepoints*trials)
    dataX = fft( dataR, nConv, 2 );
    
    % Do the fft on data-ERP (non-phaselocked part of the signal)
    % subtract the ERP from the single trials
    EEG.subdata = bsxfun(@minus,EEG.data,mean(EEG.data,3));
    subdata = double( EEG.subdata );
    subdataR = reshape(subdata, EEG.nbchan, []);
    subdataX = fft( subdataR, nConv, 2 );
       
    
    % initialize output time-frequency matrix for power and phases
    templateTF    = zeros(EEG.nbchan-2, length(frex), EEG.pnts);
    
    tmp_tf        = templateTF; %to save tf in steps in between
    tf            = templateTF; %to save final tf with power
    tf_bs         = templateTF; %to save baseline-normalized tf
    
    tmp_tf_minERP = templateTF; 
    tf_minERP     = templateTF; %to save final tf with power non-phaselocked (without ERP)
    tf_minERP_bs  = templateTF;
  
    
    % loop over frequencies
    for fi=1:num_freq
        
        % fwhm/number of cycles parameter in complex morlet wavelet:
        s = num_cycles(fi) / (2*pi*frex(fi));
        
        % create wavelet
        cmw = exp(1i * 2*pi * frex(fi) * time) .* exp( -time.^2 / (2*s^2));
        
        % FFT the wavelet kernel
        cmwX = fft( cmw, nConv );
        cmwX = cmwX./max(cmwX); %normalize the wavelet
        
        
        for chani=1:(EEG.nbchan-2)
            
            %%% the rest of convolution
            % Get the analytic signal (as)
            as = ifft( dataX(chani,:) .* cmwX ); % Multiply spectra of signal and wavelet + do inverse FFT to get to time domain
            as = as(:,(halfK+1):(end-halfK)); % cut off "wings"
            as = squeeze( reshape( as, size( EEG.data(chani,:,:) ) ) ); % cut signal based on supertrial into individual trials again
            
            % Extract result of convolution for low load + high load correct target trials:
            as = as(:,Epochs);
            
            % extract power
            aspow = mean( abs( as ).^2, 2);
            
            
            % A) Average power over trials and put power data into big matrix
            
            %%% Save mean power in timefreq matrix tf
            tmp_tf(chani,fi,:)  = aspow;            
            
            %%% For power in the non-phaselocked part of the signal (subdata):
            as2 = ifft( subdataX(chani,:) .* cmwX ); % Multiply spectra of signal and wavelet + do inverse FFT to get to time domain
            as2 = as2(:,(halfK+1):(end-halfK)); % cut off "wings"
            as2 = squeeze( reshape( as2, size( subdata(chani,:,:) ) ) ); % cut signal based on supertrial into individual trials again

            as2 = as2(:,Epochs);
            
            aspow2 = mean( abs( as2 ).^2, 2);

            %%% Save 
            tmp_tf_minERP(chani,fi,:)  = aspow2;      
            
         end
         
     end
    
    % convert baseline time into indices
    baseidx = reshape( dsearchn(EEG.times',baseline_windows'), [],2);

    % Baseline correction - db normalization of power time-freq plot
    tf_bs(:,:,:)        = 10 * log10( bsxfun(@rdivide, tmp_tf(:,:,:),        mean(tmp_tf(       :,:,baseidx(1):baseidx(2),1),3) ) );
    tf_minERP_bs(:,:,:) = 10 * log10( bsxfun(@rdivide, tmp_tf_minERP(:,:,:), mean(tmp_tf_minERP(:,:,baseidx(1):baseidx(2),1),3) ) );    
    
    % Downsample the result:
    times2save    = -500:50:2500;
    times2saveidx = dsearchn(EEG.times', times2save');
     
    tfdwnsamp(:,:,:)           = tf(          :,:,times2saveidx);
    tfdwnsamp_minERP(:,:,:)    = tf_minERP(   :,:,times2saveidx);
    tfdwnsamp_bs(:,:,:)        = tf_bs(       :,:,times2saveidx);
    tfdwnsamp_minERP_bs(:,:,:) = tf_minERP_bs(:,:,times2saveidx);
 
     
    % Save time-frequency matrix for this subject and this channel:
    savepath = '/Volumes/HER-1/onderzoeksarchief/19-452_CONTROL_BS/F_DataAnalysis/NBack EEG/Timefreq Matrix/';
    
    save([savepath 'Timefreq_power'            num2str(CurrSubj) '.mat'], 'tfdwnsamp')
    save([savepath 'Timefreq_power_minERP'     num2str(CurrSubj) '.mat'], 'tfdwnsamp_minERP')
    save([savepath 'Timefreq_power_bs_'        num2str(CurrSubj) '.mat'], 'tfdwnsamp_bs')
    save([savepath 'Timefreq_power_minERP_bs_' num2str(CurrSubj) '.mat'], 'tfdwnsamp_minERP_bs')
    
end




%% Get average values from all subjects:

% Define frontal midline channels
FMchans = [find(strcmpi( {EEG.chanlocs.labels}, 'FC1')) ...
    find(strcmpi( {EEG.chanlocs.labels}, 'FC2')) ...
    find(strcmpi( {EEG.chanlocs.labels}, 'Fz'))...
    find(strcmpi( {EEG.chanlocs.labels}, 'Cz'))];

% Define frequency of interest
numfreq = 56;
frqsvec = logspace(log10(2),log10(30), numfreq);
frqs    = [4 8]; % 4-8 Hz for theta
frqsidx = dsearchn(frqsvec', frqs');

% Define time slot of interest
tROI        = [0.2  .8]; %200-800 ms post-stim.onset
timevec     = -0.5:.05:2.5;
tROIidx     = dsearchn(timevec', tROI');
baselinewin = [ -0.35 -0.1];
bidx        = reshape( dsearchn(timevec',baselinewin'), [],2);


% initialize, theta power, seeded synchronization and IPL maps
output_thetapow = zeros(length(Subj),4);
% [synch_FM_all_lo,           synch_FM_all_hi]            = deal( zeros(length(Subj), EEG.nbchan-2, length(frqsvec), length(timevec)) );
% [FM_PLI_lo,                 FM_PLI_hi]                  = deal( zeros(length(Subj), EEG.nbchan-2, length(frqsvec), length(timevec)) );
% [FM_seeded_synch_lo,        FM_seeded_synch_hi]         = deal( zeros(length(Subj), 2,  EEG.nbchan-2) );
% [FM_seeded_synch_lo_bslnsub,FM_seeded_synch_hi_bslnsub] = deal( zeros(length(Subj), EEG.nbchan-2) );
% [FM_PLI_lo_avgROI,          FM_PLI_hi_avgROI]           = deal( zeros(length(Subj), 2,  EEG.nbchan-2) );
% [FM_PLI_lo_avgROI_bslnsub,  FM_PLI_hi_avgROI_bslnsub]   = deal( zeros(length(Subj), EEG.nbchan-2) );

% For plotting:
tf_full_all = zeros(length(Subj), EEG.nbchan-2, 56, 61);


% Change directory to output from time-frequency analysis
cd('/Volumes/HER-1/onderzoeksarchief/19-452_CONTROL_BS/F_DataAnalysis/NBack EEG/Timefreq Matrix/Non-phase-locked baseline-corrected')

% Find .mat-files in this directory
files = dir('*mat') ;
if length(Subj) ~= length(files) 
    fprintf('Missing files!')
end
N = length(files) ; 

% intialize figure windows
figure(2),clf; colormap jet
figure(3),clf; colormap jet

% Loop over subjects
for Subjecti = 1:N
    
    clear loadmat 
    
    CurrSubj = Subj(Subjecti)   
      
    %%% Time-frequency/theta power analysis 
    % Load time-freq matrix per participant:
    loadmat = load(files(Subjecti).name);
    
    % Save in outputmat
    output_thetapow(Subjecti,1) = Subj(Subjecti);
    
    % Squeeze out subject's data for Frontal midline chans only
    tf_full = loadmat.tfdwnsamp_minERP_bs;
    tf_roi  = squeeze( mean( tf_full( FMchans ,:,:), 1) );
    % Save in full matrix for average plot
    tf_full_all(Subjecti, :,:,:) = tf_full;
    
    % Create average in theta freq and time 
    tf_roi_theta = mean(tf_roi(frqsidx(1):frqsidx(2) , tROIidx(1):tROIidx(2) ),3);
    
    % Find the peak frequency for power in this theta ROI
    [maxfrqidxT, maxtimeidxT] = find( tf_roi_theta == max(max( tf_roi_theta ) ));
    maxfrqT                   = frqsvec(maxfrqidxT  + frqsidx(1) -1 );
    maxtimeT                  = timevec(maxtimeidxT + tROIidx(1) -1);
    
    % Draw box of 3 Hz x 300 ms around peak power value
    frqboxT  = dsearchn(frqsvec', [(maxfrqT-1.5)   (maxfrqT+1.5)]');
    timeboxT = dsearchn(timevec', [(maxtimeT-0.2) (maxtimeT+0.2)]');
    
    % Average power in box (ROI):
    theta = mean( mean( tf_roi( frqboxT(1):frqboxT(2) , timeboxT(1):timeboxT(2)) ) );
    
    output_thetapow(Subjecti,2) = theta;
    output_thetapow(Subjecti,3) = maxfrqT;
    output_thetapow(Subjecti,4) = maxtimeT;
    
    % Plot time-frequency
    figure(4)
    subplot(10,7,Subjecti)
    contourf(timevec,frqsvec,squeeze(mean(tf_roi(:,:,:),3)),40,'linecolor','none')
    set(gca,'clim',[0 5],'xlim',[-.2 2],'yscale','log','ytick',round(frqsvec(1:10:numfreq),0));
    xlabel('Time (s)'), ylabel('Frequency (Hz)'); c = colorbar; c.Label.String = 'Power \muV';
    rectangle('Position',[ maxtimeT-0.2  maxfrqT-1.5 0.4  3],'EdgeColor','k','linew',2)    
    
end


%%% SAVE
savepath2 = '/Volumes/HER-1/onderzoeksarchief/19-452_CONTROL_BS/F_DataAnalysis/NBack EEG/Timefreq Matrix/Non-phase-locked baseline-corrected/Outcomes/';
save([savepath2 'output_thetapow.mat'],  'output_thetapow');
writematrix(output_thetapow, [savepath2 'output_thetapow.txt']);