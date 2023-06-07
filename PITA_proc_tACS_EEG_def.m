%% Script to analyze neural oscillations in between-tACS periods

%% EEG pre-processing - PITA - all data

%% NOTE!
% Geen stimulatie uitgevoerd bij (subjectID-sessie):
% 297-2 - REAL tACS - LE, AQM   - Notes: Matlab gaf foutmelding, fenne is gekomen om te helpen, zijn ongeveer een uur later begonnen dan gepland
% 334-1 - REAL tACS - MvK, RK   - none
% 396-1 - REAL tACS - AQM, FMS  - none
% 602-1 - REAL tACS - MVK       - none
% 626-2 - REAL tACS - MvK, RK   - Notes: Voorafgaand aan tACS-EEG, verstoring signaal rondom tACS electrode door gel/vochtig (FC5/FC6/F7)
% 913-1 - REAL tACS - RK, TP    - none
% 121-2 - SHAM      - RK, MvK   - none
% 291-2 - SHAM      - RK        - none
% 362-1 - SHAM      - MvK, RK   - Notes: (post-tACS rest. EEG) Script startten niet, dus opnieuw resting-state EEG script gestart i.p.v. post-tACS EEG script. (hele experiment) Gel/wax/iets anders in haar.
% 442-1 - SHAM      - MVK, RK   - none
% 524-1 - SHAM      - AQM, RK   - none
% 559-2 - SHAM      - AQM, MvK  - Notes: Participant was vergeten het tacs script aan te zetten dus een half uur kwam ik daar achter en toen heeft ze het aangezet. [Added note: maar stimulator heeft niet gedraaid, ook niet toen script wel heeft gerund]
% 710-1 - SHAM      - Mvk, RK   - none
%
% Besluit > Excludeer ppn bij wie stimulatie niet is uitgevoerd tijdens real tACS sessie.
% Dat zijn ppn: 297, 334, 396, 602, 626, 913.
% Behoud ppn bij wie stimulatie niet is uitgevoerd in sham sessie
%
% full_subj_list = [669	557 363	638	602	989	383	502	733	442	575	913 710	262 ...
%       752 227	565	626 334	362	600	121	319 923	915	298	202	692	275	...
%       508 291	803	755	681	876	134	559	396	818	601	297	524	883	193	642];


%% Start by clearing workspace

clear
close all


%% Initialize

% open EEGlab
cd('/Users/fsmits2/Downloads/eeglab2022.1')
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

cd('/Users/fsmits2/Documents/PITA_analysis'); % return to PITA analysis folder

% set paths to data
Path2EEGsets = '/Users/fsmits2/Downloads/1 EEG data tacs-eeg/post-tacs r-s processed/';

% enter subject names
subj_list =[669	557 363	638	989	383	502	733	442	575	710	262 ...
    752 227	565	362	600	121	319 923	915	298	202	692	275	...
    508 291	803	755	681	876	134	559	818	601	524	883	193	642];
sessions  = [1 2];


% enter filenames
rec1 = 'restingstate-pretACS-';
rec2 = 'Encoding-';
rec3 = 'TACSEEG-';
rec4 = 'restingstate-posttACS-';
rec5 = 'Retrieval-';
file_type = {rec1, rec2, rec3, rec4, rec5};


%% Do the fourier Transform

% Pre-specify a matrix (dataframe) to save outcomes in the size: Subject x session x channels x condition(eyes open vs. closed) x outcome(theta power, beta power, theta/beta ratio, number of 1-sec-epochs)
% % % datfr       = nan(length(subj_list), length(sessions), 30, 513, 600); %for tACS-EEG data
% % % datfr_dense = datfr;
datfr_rso   = nan(length(subj_list), length(sessions), 30, 513, 120); %for resting-state EEG data eyes-open
datfr_rsc   = nan(length(subj_list), length(sessions), 30, 513, 120); %eyes-closed

% Which task (file type) you want to analyze?
fileno = 4;

% Loop over subjects and sessions
for subj_i = 15:length(subj_list)
    for sess_i = 1:length(sessions)

        fprintf('\n*** Load preprocessed data from: subject %i session %i\n', subj_list(subj_i), sessions(sess_i));
        fileName = [file_type{fileno} num2str(subj_list(subj_i)) '-' num2str(sess_i) '_CleanEEG.set'];

        % Load EEG set
        EEG      = pop_loadset('filename', fileName, 'filepath', Path2EEGsets);

        % Split dataset into eyes-open and eyes-closed
        openidx   = find(11 == str2double({EEG.event.type}));
        closedidx = find(22 == str2double({EEG.event.type}));
        openeps   = cell2mat({EEG.event(openidx).epoch});
        closedeps = cell2mat({EEG.event(closedidx).epoch});
        % Put EEG data in double precision for good computation performance
% % %         EEG.datax = double( EEG.data); 
        EEG.datax_o = double( EEG.data(:,:,openeps) ); 
        EEG.datax_c = double( EEG.data(:,:,closedeps) ); 

        % Prepare Fourier Transform
        nfft      = EEG.srate * 4; % For zero-padding (upsampling) and overlapping in Welch's method
        nOverlap  = size(EEG.datax_o,2)/2; % For no overlap: 0;  for 50% overlap: size(EEG.datax,2)/2

        hannw     = hann(EEG.pnts); % Create Hann window to taper the data with.      = .5 * (1 - cos(2*pi*linspace( 0, 1, size(EEG.datax,2) ) ));

        powspec_o   = nan( EEG.nbchan-2, EEG.srate*2+1, size(EEG.datax_o,3) ); % pre-specify power spectrum variable to save FFT results
        powspec_c   = nan( EEG.nbchan-2, EEG.srate*2+1, size(EEG.datax_c,3) ); % pre-specify power spectrum variable to save FFT results

        % Do the FFT. Loop over channels
        fprintf('\n*** Compute power spectrum - subject %i session %i\n', subj_list(subj_i), sessions(sess_i));
        for chani = 1:EEG.nbchan-2 % minus the last two channels (EMG): HEOG & VEOG
            % Do the FFT for each frequency and epoch using Welch's method (MATLAB's function pwelch)
            [powspec_o(chani,:,:), hz] = pwelch( squeeze( EEG.datax_o(chani,:,:) ), hannw, nOverlap, nfft, EEG.srate );
            [powspec_c(chani,:,:), hz] = pwelch( squeeze( EEG.datax_c(chani,:,:) ), hannw, nOverlap, nfft, EEG.srate );
        end

        % Save powespectra in the dense matrix:
        fprintf('\n*** Save dense power spectrum - subject %i session %i\n', subj_list(subj_i), sessions(sess_i));
        datfr_rso(subj_i, sess_i, :, : , 1:size(EEG.datax_o,3)) = powspec_o;
        datfr_rsc(subj_i, sess_i, :, : , 1:size(EEG.datax_c,3)) = powspec_c;

% % %         %%% Only for tACS-EEG
% % %         % Save powerspectra per epoch in it's original position in the tACS-EEG timing
% % %         fprintf('\n*** Save power spectrum per epoch in original timing - subject %i session %i\n', subj_list(subj_i), sessions(sess_i));
% % %         for epoci = 1:EEG.trials      
% % %             timetrigg = str2double( EEG.event(epoci).type ) ;
% % %             block     = floor( timetrigg );
% % %             sec       = round( (timetrigg - block) * 100 );
% % % 
% % %             epocpos   = (block-1) * 29 + sec;
% % %             datfr(subj_i, sess_i, :, : , epocpos) = powspec(:, :, epoci);
% % %         end

        clear EEG
        ALLEEG(1:end) = [];
        clear powspec_o; clear EEG.datax_o; clear powspec_c; clear EEG.datax_c;

    end
end

% Write to file
filename2='hz_saved.mat';
save(filename2,'hz');

filename='datfr_rso_saved.mat';
save(filename,'datfr_rso', '-v7.3');

filename1='datfr_rsc_saved.mat';
save(filename1,'datfr_rsc', '-v7.3');

% % % filename='datfr_dense_saved.mat';
% % % save(filename,'datfr_dense', '-v7.3');



%% Look at result | Topoplots for channel selection in spectral power analyses

% Load the table with noisy channels (to be excluded from analysis)
bdchns = table2cell(  readtable( [Path2EEGsets '/Overview_badchannels_05-Jun-2023.txt'] ,'Format','auto') );

% pre-specify variable to save clean power outcomes in
% % % cleandatfr    = nan(length(subj_list), length(sessions), 30, 513, 600); % for tacs-EEG data
cleandatfr_o    = nan(length(subj_list), length(sessions), 30, 513, 120); %eyes-open
cleandatfr_c    = nan(length(subj_list), length(sessions), 30, 513, 120); %eyes-closed

% Load EEG set
fileName = [file_type{1} '642-2_CleanEEG.set'];
EEG      = pop_loadset('filename', fileName, 'filepath', Path2EEGsets);

for subj_i = 1:length(subj_list)
    for sess_i = 1:length(sessions)

        chanarray = 1:30;

        % Select only from channels with clean data, i.e. remove noisy channels from analysis outcome:
        badchans     = regexp(bdchns{subj_i,sess_i+1} , ',', 'split');
        badchanarray = [];
        if sum( strcmpi( badchans, '0') ) < 1
            for bchni = 1:length(badchans)
                badchanarray(bchni) = find( strcmpi( badchans{bchni}, {EEG.chanlocs.labels} ));
            end
        end
        chanarray(badchanarray) = []; % remove the noisy channels from chanarray

        % Copy only from clean channels data from alldatfr to powdatfr
        fprintf('\n*** Save powspec from only no-noise channels from: subject %i session %i\n', subj_list(subj_i), sessions(sess_i));
        cleandatfr_o(subj_i, sess_i, chanarray, :, :) = datfr_rso(subj_i, sess_i, chanarray, :, :);
        cleandatfr_c(subj_i, sess_i, chanarray, :, :) = datfr_rsc(subj_i, sess_i, chanarray, :, :);

    end
end

% filename='cleandatfr_saved.mat';
% save(filename,'cleandatfr', '-v7.3');

filename='cleandatfr_o_saved.mat';
save(filename,'cleandatfr_o', '-v7.3');

filename1='cleandatfr_c_saved.mat';
save(filename1,'cleandatfr_c', '-v7.3');


% Find Individual Alpha Peak Frequency (IAF)
maxfrq_c = nan(length(subj_list), length(sessions), 1);
for subj_i = 1:length(subj_list)
    for sess_i = 1:length(sessions)

%        datA = log10( squeeze( mean(mean(cleandatfr_o(subj_i,sess_i,13:17,:,:),5,'omitnan'),3,'omitnan') ) );
        datB = log10( squeeze( mean(mean(cleandatfr_c(subj_i,sess_i,13:17,:,:),5,'omitnan'),3,'omitnan') ) );
%         figure(3)
%         plot(hz, datA,'k','LineWidth',3); hold on
%         plot(hz, datB,'r','LineWidth',3)
        alphafrqidx   = dsearchn(hz, [8  13]');
        frqsvec = alphafrqidx(1):alphafrqidx(2);
        datB2 = datB(alphafrqidx(1):alphafrqidx(2));
        [maxpow, idx]  = max( datB2 );
        maxfrqidx      = frqsvec(idx(1));
        maxfrq         = hz(maxfrqidx);

        maxfrq_c(subj_i, sess_i, :) = maxfrq;

    end
end

filename3='maxfrq_c_saved.mat';
save(filename3,'maxfrq_c');


% Select frequencies
frqidx   = dsearchn(hz, [4  7]');

% Average over selected channels and frequencies
powplot = nan(length(subj_list), length(sessions));

for subj_i = 1:length(subj_list)
    for sess_i = 1:length(sessions)

        % Remove noisy
        chanarray    = 1:EEG.nbchan-2;
        badchans     = regexp(bdchns{subj_i,sess_i+1} , ',', 'split');
        badchanarray = [];
        if sum( strcmpi( badchans, '0') ) < 1
            for bchni = 1:length(badchans)
                badchanarray(bchni) = find( strcmpi( badchans{bchni}, {EEG.chanlocs.labels} ));
            end
        end
        chanarray(badchanarray) = []; % remove the noisy channels from chanarray
        
        % Extract from data
        dataextract = [];
        dataextract = datfr(subj_i,sess_i, chanarray, frqidx(1):frqidx(2), :);
        % Average over channels and frequencies and log-transform
        powplot(subj_i,sess_i) = log10(squeeze( mean( mean(dataextract,5,'omitnan') ,4,'omitnan') ) );
       
    end
end

% Define frequency bands
betawin  = [14  30];
betaidx  = dsearchn(hz, betawin');
alphawin = [8  13];
alphaidx = dsearchn(hz, alphawin');
thetawin = [4  7.5];
thetaidx = dsearchn(hz, thetawin');
frqidx   = dsearchn(hz, [4.5  5.5]');

% Extract average power in each band over all epochs, subjects and sessions
powtopobeta  = log10( squeeze( mean( mean( mean( mean( cleandatfr(:, :, :, betaidx(1):betaidx(2),   :) ,5,'omitnan') ,4,'omitnan') ,2,'omitnan') ,1,'omitnan') ) );
powtopoalpha = log10( squeeze( mean( mean( mean( mean( cleandatfr(:, :, :, alphaidx(1):alphaidx(2), :) ,5,'omitnan') ,4,'omitnan') ,2,'omitnan') ,1,'omitnan') ) );
powtopotheta = log10( squeeze( mean( mean( mean( mean( cleandatfr(:, :, :, thetaidx(1):thetaidx(2), :) ,5,'omitnan') ,4,'omitnan') ,2,'omitnan') ,1,'omitnan') ) );
powtopo5hz   = log10( squeeze( mean( mean( mean( mean( cleandatfr(:, :, :, frqidx(1):frqidx(2),     :) ,5,'omitnan') ,4,'omitnan') ,2,'omitnan') ,1,'omitnan') ) );

% Plot
figure(2); colormap jet
subplot(1,4,1); plot( hz, powspec )
subplot(1,4,2); topoplotIndie(powtopobeta,  EEG.chanlocs(1:30),'numcontour',0,'electrodes','numbers','shading','interp');
title({'14 - 30 Hz'}); c = colorbar; c.Label.String = 'Power \muV^2'; set(gca,'clim',[-0.3 0.2]);
subplot(1,4,3); topoplotIndie(powtopoalpha, EEG.chanlocs(1:30),'numcontour',0,'electrodes','numbers','shading','interp');
title({'8 - 13 Hz'});    c = colorbar; c.Label.String = 'Power \muV^2'; set(gca,'clim',[0.2 0.9]);
subplot(1,4,4); topoplotIndie(powtopotheta, EEG.chanlocs(1:30),'numcontour',0,'electrodes','numbers','shading','interp');
title({'4 - 7.5 Hz'});   c = colorbar; c.Label.String = 'Power \muV^2'; set(gca,'clim',[0.2 0.9]);


%% Output matrix for stats in R

% LOAD datfr as matrix (not cleandatfr)
load hz_saved.mat
% % % datfr = importdata('datfr_saved.mat');
datfr_rso = importdata('datfr_rso_saved.mat');
datfr_rsc = importdata('datfr_rsc_saved.mat');

% Load the table with noisy channels (to be excluded from analysis)
bdchns = table2cell(  readtable( [Path2EEGsets '/Overview_badchannels_14-Feb-2023.txt'] ,'Format','auto') );

% Load EEG set
fileName = [file_type{1} '642-2_CleanEEG.set'];
EEG      = pop_loadset('filename', fileName, 'filepath', Path2EEGsets);

% Select channels for. For theta activity: overall power is strongest over Fz and surrounding channels
channames = {'AF3' 'AF4' 'FC1' 'FC2' 'FC5' 'FC6' 'F7' 'F8'};
chans2use = [];
for chani = 1:length(channames)
    chans2use(chani) = find( strcmpi( channames{chani}, {EEG.chanlocs.labels} ));
end

% Average over selected channels and frequencies
% % % powdat = nan(length(subj_list), length(sessions), 600); %for tACS_EEG data
powdat = nan(length(subj_list), length(sessions), 120); %eyes-open

% Select frequencies
frqidx   = dsearchn(hz, [4.75  5.25]');

for subj_i = 11:length(subj_list)
    for sess_i = 1:length(sessions)

        % Remove noisy
        chanarray    = 1:EEG.nbchan-2;
        badchans     = regexp(bdchns{subj_i,sess_i+1} , ',', 'split');
        badchanarray = [];
        if sum( strcmpi( badchans, '0') ) < 1
            for bchni = 1:length(badchans)
                badchanarray(bchni) = find( strcmpi( badchans{bchni}, {EEG.chanlocs.labels} ));
            end
        end
        chanarray(badchanarray) = []; % remove the noisy channels from chanarray
        chans2use_now = chans2use(ismember(chans2use,chanarray));
        
        % Extract from data
        dataextract = [];
        dataextract = datfr_rso(subj_i,sess_i, chans2use_now, frqidx(1):frqidx(2), :);
        % Average over channels and frequencies and log-transform
        powdat(subj_i,sess_i,:) = log10(squeeze( mean( mean(dataextract,4,'omitnan') ,3,'omitnan') ) );
       
    end
end

filename='powdat_rso_5Hz_saved.mat';
save(filename,'powdat', '-v6');



%% Plots for BRST2023 poster

% LOAD datfr as matrix
load hz_saved.mat
datfr = importdata('datfr_saved.mat');

% remove data from datasets with gelbridge - subject 989 (session 1) and 818 (session 2)
subj989 = find(subj_list==989); datfr(subj989,1,:,:,:) = NaN;
subj818 = find(subj_list==818); datfr(subj818,2,:,:,:) = NaN;

% Load the table with noisy channels (to be excluded from analysis)
bdchns = table2cell(  readtable( [Path2EEGsets '/Overview_badchannels_20-Jan-2023.txt'] ,'Format','auto') );

% Load EEG set
fileName = 'TACSEEG-669-1_CleanEEG_CleanEpochs.set';
EEG      = pop_loadset('filename', fileName, 'filepath', Path2EEGsets);

% Average over selected channels and frequencies
powdatfull = nan(length(subj_list), length(sessions), 30, 600);

% Select frequencies
frqidx   = dsearchn(hz, [4.75  5.25]');

for subj_i = 1:length(subj_list)
    for sess_i = 1:length(sessions)

        % Remove noisy channels
        chanarray    = 1:EEG.nbchan-2;
        badchans     = regexp(bdchns{subj_i,sess_i+1} , ',', 'split');
        badchanarray = [];
        if sum( strcmpi( badchans, '0') ) < 1
            for bchni = 1:length(badchans)
                badchanarray(bchni) = find( strcmpi( badchans{bchni}, {EEG.chanlocs.labels} ));
            end
        end
        chanarray(badchanarray) = [];

        % Extract from data per (clean) channel
        for chani = chanarray
            dataextractfull = [];
            dataextractfull = datfr(subj_i, sess_i, chani, frqidx(1):frqidx(2), :);
            % Average over channels and frequencies and log-transform
            powdatfull(subj_i,sess_i,chani,:) = squeeze( mean(dataextractfull,4,'omitnan') );
        end

    end
end

% Read in conditions (real vs. sham tACS)
conditions = readtable('Random_allocation_log_PITA_deblinded_matlab.csv');
subjects = table2array(conditions(:,1));
subjmatches = find( ismember(subjects, subj_list) > 0 );
conds = conditions(subjmatches,:);
% Couple conditions to dataframe rows
s1_real = find( table2array(conds(:,2)) >  0 )';
s1_sham = find( table2array(conds(:,2)) == 0 )';
s2_real = find( table2array(conds(:,3)) >  0 )';
s2_sham = find( table2array(conds(:,3)) == 0 )';

% Define where to cut trials to divide data into 5 segments (of 4x 30-sec EEG recording each)
cuts = 1:116:580;

% Data per segment per condition
powtopotheta1_sham  = mean( [ mean( log10( squeeze( mean( powdatfull(s1_sham, 1, :, cuts(1):cuts(2)-1) ,4,'omitnan') ) ) ,1,'omitnan') ...
                        mean( log10( squeeze( mean( powdatfull(s2_sham, 2, :, cuts(1):cuts(2)-1) ,4,'omitnan') ) ) ,1,'omitnan') ] ,1) ;
powtopotheta2_sham  = mean( [ mean( log10( squeeze( mean( powdatfull(s1_sham, 1, :, cuts(2):cuts(3)-1) ,4,'omitnan') ) ) ,1,'omitnan') ...
                         mean( log10( squeeze( mean( powdatfull(s2_sham, 2, :, cuts(2):cuts(3)-1) ,4,'omitnan') ) ),1,'omitnan') ] ,1) ;
powtopotheta3_sham  = mean( [ mean( log10( squeeze( mean( powdatfull(s1_sham, 1, :, cuts(3):cuts(4)-1) ,4,'omitnan') ) ) ,1,'omitnan') ...
                        mean( log10( squeeze( mean( powdatfull(s2_sham, 2, :, cuts(3):cuts(4)-1) ,4,'omitnan') ) ),1,'omitnan') ] ,1) ;
powtopotheta4_sham  = mean( [ mean( log10( squeeze( mean( powdatfull(s1_sham, 1, :, cuts(4):cuts(5)-1) ,4,'omitnan') ) ) ,1,'omitnan') ...
                        mean( log10( squeeze( mean( powdatfull(s2_sham, 2, :, cuts(4):cuts(5)-1) ,4,'omitnan') ) ),1,'omitnan') ] ,1) ;
powtopotheta5_sham  = mean( [ mean( log10( squeeze( mean( powdatfull(s1_sham, 1, :, cuts(5):580)       ,4,'omitnan') ) ) ,1,'omitnan') ...
                        mean( log10( squeeze( mean( powdatfull(s2_sham, 2, :, cuts(5):580)       ,4,'omitnan') ) ),1,'omitnan') ] ,1) ;
powtopotheta1_real  = mean( [ mean( log10( squeeze( mean( powdatfull(s1_real, 1, :, cuts(1):cuts(2)-1) ,4,'omitnan') ) ) ,1,'omitnan') ...
                        mean( log10( squeeze( mean( powdatfull(s2_real, 2, :, cuts(1):cuts(2)-1) ,4,'omitnan') ) ),1,'omitnan') ] ,1) ;
powtopotheta2_real  = mean( [ mean( log10( squeeze( mean( powdatfull(s1_real, 1, :, cuts(2):cuts(3)-1) ,4,'omitnan') ) ) ,1,'omitnan') ...
                        mean( log10( squeeze( mean( powdatfull(s2_real, 2, :, cuts(2):cuts(3)-1) ,4,'omitnan') ) ),1,'omitnan') ] ,1) ;
powtopotheta3_real  = mean( [ mean( log10( squeeze( mean( powdatfull(s1_real, 1, :, cuts(3):cuts(4)-1) ,4,'omitnan') ) ) ,1,'omitnan') ...
                        mean( log10( squeeze( mean( powdatfull(s2_real, 2, :, cuts(3):cuts(4)-1) ,4,'omitnan') ) ),1,'omitnan') ] ,1) ;
powtopotheta4_real  = mean( [ mean( log10( squeeze( mean( powdatfull(s1_real, 1, :, cuts(4):cuts(5)-1) ,4,'omitnan') ) ) ,1,'omitnan') ...
                        mean( log10( squeeze( mean( powdatfull(s2_real, 2, :, cuts(4):cuts(5)-1) ,4,'omitnan') ) ),1,'omitnan') ] ,1) ;
powtopotheta5_real  = mean( [ mean( log10( squeeze( mean( powdatfull(s1_real, 1, :, cuts(5):580)       ,4,'omitnan') ) ) ,1,'omitnan') ...
                        mean( log10( squeeze( mean( powdatfull(s2_real, 2, :, cuts(5):580)       ,4,'omitnan') ) ),1,'omitnan') ] ,1) ;

% Data per segment per condition
powtopotheta1_sham  = mean( log10( squeeze( mean( powdatfull(s2_sham, 2, :, cuts(1):cuts(2)-1) ,4,'omitnan') ) ),1,'omitnan');
powtopotheta2_sham  = mean( log10( squeeze( mean( powdatfull(s2_sham, 2, :, cuts(2):cuts(3)-1) ,4,'omitnan') ) ),1,'omitnan') ;
powtopotheta3_sham  = mean( log10( squeeze( mean( powdatfull(s2_sham, 2, :, cuts(3):cuts(4)-1) ,4,'omitnan') ) ),1,'omitnan') ;
powtopotheta4_sham  = mean( log10( squeeze( mean( powdatfull(s2_sham, 2, :, cuts(4):cuts(5)-1) ,4,'omitnan') ) ),1,'omitnan') ;
powtopotheta5_sham  = mean( log10( squeeze( mean( powdatfull(s2_sham, 2, :, cuts(5):580)       ,4,'omitnan') ) ),1,'omitnan') ;
powtopotheta1_real  = mean( log10( squeeze( mean( powdatfull(s2_real, 2, :, cuts(1):cuts(2)-1) ,4,'omitnan') ) ),1,'omitnan') ;
powtopotheta2_real  = mean( log10( squeeze( mean( powdatfull(s2_real, 2, :, cuts(2):cuts(3)-1) ,4,'omitnan') ) ),1,'omitnan') ;
powtopotheta3_real  = mean( log10( squeeze( mean( powdatfull(s2_real, 2, :, cuts(3):cuts(4)-1) ,4,'omitnan') ) ),1,'omitnan') ;
powtopotheta4_real  = mean( log10( squeeze( mean( powdatfull(s2_real, 2, :, cuts(4):cuts(5)-1) ,4,'omitnan') ) ),1,'omitnan') ;
powtopotheta5_real  = mean( log10( squeeze( mean( powdatfull(s2_real, 2, :, cuts(5):580)       ,4,'omitnan') ) ),1,'omitnan') ;



figure(6); colormap jet
clims_theta = [0.2 0.9];
subplot(6,2,1); topoplot( mean( [powtopotheta1_real; powtopotheta2_real; powtopotheta3_real; powtopotheta4_real; powtopotheta5_real ],1), EEG.chanlocs(1:30)); %,'numcontour',0,'electrodes','numbers','shading','interp');
title({'Real tACS, overall'}); c = colorbar; c.Label.String = 'Power \muV^2 (log10)'; set(gca,'clim',clims_theta);
subplot(6,2,3); topoplot(powtopotheta1_real,  EEG.chanlocs(1:30)); %,'numcontour',0,'electrodes','numbers','shading','interp');
%title({'1-4 min. real tACS'}); 
set(gca,'clim',clims_theta);
subplot(6,2,5); topoplot(powtopotheta2_real, EEG.chanlocs(1:30)); %,'numcontour',0,'electrodes','numbers','shading','interp');
%title({'5-8 min. real tACS'}); 
set(gca,'clim',clims_theta);
subplot(6,2,7); topoplot(powtopotheta3_real, EEG.chanlocs(1:30)); %,'numcontour',0,'electrodes','numbers','shading','interp');
%title({'9-12 min. real tACS'}); 
set(gca,'clim',clims_theta);
subplot(6,2,9); topoplot(powtopotheta4_real, EEG.chanlocs(1:30)); %,'numcontour',0,'electrodes','numbers','shading','interp');
%title({'13-16 min. real tACS'}); 
set(gca,'clim',clims_theta);
subplot(6,2,11); topoplot(powtopotheta5_real, EEG.chanlocs(1:30)); %,'numcontour',0,'electrodes','numbers','shading','interp');
%title({'17-20 min. real tACS'}); 
set(gca,'clim',clims_theta);
subplot(6,2,2); topoplot( mean( [powtopotheta1_sham; powtopotheta2_sham; powtopotheta3_sham; powtopotheta4_sham; powtopotheta5_sham ],1), EEG.chanlocs(1:30)); %,'numcontour',0,'electrodes','numbers','shading','interp');
title({'Sham, overall'}); c = colorbar; c.Label.String = 'Power \muV^2 (log10)'; set(gca,'clim',clims_theta);
subplot(6,2,4); topoplot(powtopotheta1_sham,  EEG.chanlocs(1:30)); %,'numcontour',0,'electrodes','numbers','shading','interp');
%title({'1-4 min. sham'}); 
set(gca,'clim',clims_theta);
subplot(6,2,6); topoplot(powtopotheta2_sham, EEG.chanlocs(1:30)); %,'numcontour',0,'electrodes','numbers','shading','interp');
%title({'5-8 min. sham'}); 
set(gca,'clim',clims_theta);
subplot(6,2,8); topoplot(powtopotheta3_sham, EEG.chanlocs(1:30)); %,'numcontour',0,'electrodes','numbers','shading','interp');
%title({'9-12 min. sham'}); 
set(gca,'clim',clims_theta);
subplot(6,2,10); topoplot(powtopotheta4_sham, EEG.chanlocs(1:30)); %,'numcontour',0,'electrodes','numbers','shading','interp');
%title({'13-16 min. sham'}); 
set(gca,'clim',clims_theta);
subplot(6,2,12); topoplot(powtopotheta5_sham, EEG.chanlocs(1:30)); %,'numcontour',0,'electrodes','numbers','shading','interp');
%title({'17-20 min. sham'}); 
set(gca,'clim',clims_theta);

set(gcf, 'color', 'none');    
set(gca, 'color', 'none');



% plot real-sham difference maps
diffmaptheta1_sham = powtopotheta1_real - powtopotheta1_sham;
diffmaptheta2_sham = powtopotheta2_real - powtopotheta2_sham;
diffmaptheta3_sham = powtopotheta3_real - powtopotheta3_sham;
diffmaptheta4_sham = powtopotheta4_real - powtopotheta4_sham;
diffmaptheta5_sham = powtopotheta5_real - powtopotheta5_sham;

% Select channels for. For theta activity: overall power is strongest over Fz and surrounding channels
channames = {'AF3' 'AF4' 'FC1' 'FC2' 'FC5' 'FC6' 'F7' 'F8'};
chans2use = [];
for chani = 1:length(channames)
    chans2use(chani) = find( strcmpi( channames{chani}, {EEG.chanlocs.labels} ));
end

figure(8); colormap jet
clims_diff_theta = [-0.15 0.15];
subplot(1,5,1); topoplot(diffmaptheta1_sham,  EEG.chanlocs(chans2use));
title({'real-sham difference, theta power (4-7 Hz), s1'});    c = colorbar; c.Label.String = 'Power \muV^2'; set(gca,'clim',clims_diff_theta);
subplot(1,5,2); topoplot(diffmaptheta2_sham,   EEG.chanlocs(chans2use));
title({'real-sham difference, theta power (4-7 Hz), s2'});    c = colorbar; c.Label.String = 'Power \muV^2'; set(gca,'clim',clims_diff_theta);
subplot(1,5,3); topoplot(diffmaptheta3_sham,   EEG.chanlocs(chans2use));
title({'real-sham difference, theta power (4-7 Hz), s3'});    c = colorbar; c.Label.String = 'Power \muV^2'; set(gca,'clim',clims_diff_theta);
subplot(1,5,4); topoplot(diffmaptheta4_sham,   EEG.chanlocs(chans2use));
title({'real-sham difference, theta power (4-7 Hz), s4'});    c = colorbar; c.Label.String = 'Power \muV^2'; set(gca,'clim',clims_diff_theta);
subplot(1,5,5); topoplot(diffmaptheta5_sham,   EEG.chanlocs(chans2use));
title({'real-sham difference, theta power (4-7 Hz), s5'});    c = colorbar; c.Label.String = 'Power \muV^2'; set(gca,'clim',clims_diff_theta);








%% Cross-frequency coupling

% Load the table with noisy channels (to be excluded from analysis)
bdchns = table2cell(  readtable( [Path2EEGsets '/Overview_badchannels_20-Jan-2023.txt'] ,'Format','auto') );

% Load EEG set
fileName = 'TACSEEG-642-2_CleanEEG_CleanEpochs.set';
EEG      = pop_loadset('filename', fileName, 'filepath', Path2EEGsets);

% Select channels for. For theta activity: overall power is strongest over Fz and surrounding channels
channames = {'AF3' 'AF4' 'FC1' 'FC2' 'FC5' 'FC6' 'F7' 'F8'};
chans2use = [];
for chani = 1:length(channames)
    chans2use(chani) = find( strcmpi( channames{chani}, {EEG.chanlocs.labels} ));
end

% Sensors to combine
combinedSensors  = chans2use;
NcombinedSensors = length(combinedSensors);
% Sampling rate [Hz]
f = EEG.srate;
% Sampling time [s]
T = 1/f;
% The number of samples to cut off from the start and end of every epoch. Ideally equal to the filter order of the lower frequency.
edgeLength = 18;
% Epoch length [s]
EpochLength = EEG.xmax-EEG.xmin;
% Number of samples per epoch
epochSamples    = EpochLength * f;
epochSamplesCut = epochSamples - 2*edgeLength;
% Time vector [s]
t = 0:T:EpochLength-T;
% Number of sensors
sensors = EEG.nbchan;
% Number of permutations
permutations = 500;
% Frequency band range [Hz]
betaLimits  = [14, 30];
thetaLimits = [4,  7];

export          = zeros(length(subj_list),2,NcombinedSensors);
export_average  = zeros(length(subj_list),2,1);
exportZ         = export;
exportZ_average = zeros(length(subj_list),2,1);
exportthetapow  = export;
exportbetapow   = export;

% DIT MOET NOG ERGENS ---------
%  % Remove noisy channels
%         chanarray    = 1:EEG.nbchan-2;
%         badchans     = regexp(bdchns{subj_i,sess_i+1} , ',', 'split');
%         badchanarray = [];
%         if sum( strcmpi( badchans, '0') ) < 1
%             for bchni = 1:length(badchans)
%                 badchanarray(bchni) = find( strcmpi( badchans{bchni}, {EEG.chanlocs.labels} ));
%             end
%         end
%         chanarray(badchanarray) = []; % remove the noisy channels from chanarray     
%       chans2use = chans2use(ismember(chans2use,chanarray));


for subj_i = 1:length(subj_list)
    for sess_i = 1:length(sessions)

        fprintf('\n*** Load tACS-EEG preprocessed data from: subject %i session %i\n', subj_list(subj_i), sessions(sess_i));
        fileName = [file_type{3} num2str(subj_list(subj_i)) '-' num2str(sess_i) '_CleanEEG_CleanEpochs.set'];

        % Load EEG set
        EEG      = pop_loadset('filename', fileName, 'filepath', Path2EEGsets);
   
        % Number of epochs
        epochs         = EEG.trials;
        % Number of samples per file
        fileSamples    = epochSamples * epochs;
        fileSamplesCut = epochSamplesCut * epochs;

        % Phase variable initialisation
        phaseCombined  = zeros(fileSamplesCut*NcombinedSensors,1);

        % Theta and beta band amplitude variable initialisation
        amplitudeThetaCombined = zeros(fileSamplesCut*NcombinedSensors,1);
        amplitudeCombined      = zeros(fileSamplesCut*NcombinedSensors,1);

        % Cycle through sensors
        for sensorIndex = 1:NcombinedSensors

            fprintf('Subject %i, Session %i, Sensor %d / %d\n',subj_list(subj_i),sess_i,sensorIndex,NcombinedSensors)

            % FMS - IIR Butterworth filter bandpass on thetalimits, filter order 4 (-48 dB rolloff).
            theta  = pop_basicfilter( EEG, combinedSensors, 'Cutoff',  thetaLimits, 'Design', 'butter', 'Filter', 'bandpass', 'Order',  8 );
            beta   = pop_basicfilter( EEG, combinedSensors , 'Cutoff',  betaLimits, 'Design', 'butter', 'Filter', 'bandpass', 'Order',  8);
            % Ik denk dat je in plaats hiervan gewoon de respectievelijk
            % theta indices en beta indices uit moet halen, en daarop je
            % Hilbert doen. Kijk Hilbert nog een x van Mike. Nu snap ik nl.
            % ook niet waarom je die start en end dingen nodig hebt. Ze
            % zitten in het 'domein' van de frequenties - of wordt dat
            % omgezet?


            % Theta band phase variable initialisation
            phase          = zeros(fileSamplesCut,1);
            amplitudeTheta = zeros(fileSamplesCut,1);

            % Beta band amplitude variable initialisation
            amplitude = zeros(fileSamplesCut,1);

            % Cycle through the epochs
            % Cycle through 12 random epochs
            ranepochs =  randi(epochs,1,12);

            for epochi = ranepochs
                % Start and end indices for the current original epoch
                iStart = (epochi-1)*epochSamples + 1;
                iEnd   = epochi*epochSamples;

                % Start and end indices for the current cut epoch
                iStartCut = (epochi-1)*epochSamplesCut + 1;
                iEndCut   = epochi*epochSamplesCut;

                % Theta band phase + amplitude and beta amplitude through the Hilbert transform
                phaseEpoch          = angle(hilbert(theta.data(sensorIndex,:,epochi)));
                amplitudeThetaEpoch = abs(hilbert(  theta.data(sensorIndex,:,epochi)));
                amplitudeEpoch      = abs(hilbert(  beta.data( sensorIndex,:,epochi)));

                % Cut the edge effects from the Hilbert-transformed epoch.
                phase(iStartCut:iEndCut)          = phaseEpoch(         edgeLength + 1:end - edgeLength);
                amplitudeTheta(iStartCut:iEndCut) = amplitudeThetaEpoch(edgeLength + 1:end - edgeLength);
                amplitude(iStartCut:iEndCut)      = amplitudeEpoch(     edgeLength + 1:end - edgeLength);
            end

% 
%             % dPAC
%             % Start and end indices for current sensor
%             iStart = (sensorIndex-1)*fileSamplesCut + 1;
%             iEnd   = sensorIndex*fileSamplesCut;
% 
%             % Append sensors data, to be able to calculate average
%             phaseCombined(iStart:iEnd)          = phase;
%             amplitudeCombined(iStart:iEnd)      = amplitude;
%             amplitudeThetaCombined(iStart:iEnd) = amplitudeTheta;
% 
%             % Phase Clustering bias
%             PCbias = mean(exp(1i*phase));
% 
%             % Debiased Phase-Ampltiude Cross-Frequency Coupling (dPAC)
%             dPAC = mean((exp(1i*phase) - PCbias) .* amplitude);
% 
%             % Save result in matrix
%             export(subj_i,sess_i,sensorIndex) = abs(dPAC);
% 
%             % Null for z-value
%             dPACnull = zeros(1,permutations);
% 
%             % Permutate the signal
%             for permutation = 1:permutations
% 
%                 % -- cut-and-paste a random portion of the data; this preserves temporal autocorrelation while removing the coupling
%                 cutLoc = 5 + randperm(length(phase)-10); % -- 5 and 10 prevent the first and last time points from being selected
%                 cutLoc = cutLoc(1);
%                 phaseShuffled = phase([cutLoc:end 1:cutLoc-1]);
% 
%                 % Compute surrogate dPAC
%                 dPACnull(permutation) = abs(mean((exp(1i*phaseShuffled) - mean(exp(1i*phaseShuffled))) .* amplitude));
%             end
% 
%             % dPAC z-score
%         meandPACnull = mean(dPACnull);
%         stddPACnull  = std(dPACnull);
%         dPACz        = (abs(dPAC) - meandPACnull) ./ stddPACnull;
% 
%         % Save result in array
%         exportZ(subj_i,sess_i,sensorIndex) = abs(dPACz);
% 
%         % FMS - Save amplitude results in array:
%         exportthetapow(subj_i,sess_i,sensorIndex) = mean(amplitudeThetaCombined)^2;
%         exportbetapow(subj_i,sess_i,sensorIndex)  = mean(amplitudeCombined)^2;
%         
        end

        clear EEG
        ALLEEG(1:end) = [];

    end
end

% dPAC for the average of the sensors
% Average the dPAC and Z-value per participant
export_average = mean(export,2);
exportZ_average = mean(exportZ,2);



%% Plot

% plot
pows         = log10(       mean( powspec(chani,        :   , : ) ,3)    );
powtop_gamma = log10( mean( mean( powspec(chanarray, 141:213, : ) ,3),2) );
powtop_alpha = log10( mean( mean( powspec(chanarray,  33:53,  : ) ,3),2) );
figure(2);  colormap jet
subplot(2,3,1); plot( hz(5:140), pows(1,5:140) )
subplot(2,3,2); topoplotIndie(powtop_gamma, EEG.chanlocs(chanarray),'numcontour',0,'electrodes','numbers','shading','interp');
subplot(2,3,3); topoplotIndie(powtop_alpha, EEG.chanlocs(chanarray),'numcontour',0,'electrodes','numbers','shading','interp');

pows_LP         = log10(       mean( powspec_LP(chani,        :   , : ) ,3)    );
powtop_gamma_LP = log10( mean( mean( powspec_LP(chanarray, 141:213, : ) ,3),2) );
powtop_alpha_LP = log10( mean( mean( powspec_LP(chanarray,  33:53,  : ) ,3),2) );
subplot(2,3,4); plot( hz(5:140), pows_LP(1,5:140), 'k-' )
subplot(2,3,5); topoplotIndie(powtop_gamma_LP, EEG.chanlocs(chanarray),'numcontour',0,'electrodes','numbers','shading','interp');
subplot(2,3,6); topoplotIndie(powtop_alpha_LP, EEG.chanlocs(chanarray),'numcontour',0,'electrodes','numbers','shading','interp');


% plot 2
channames = {'AF3' 'AF4' 'FC1' 'FC2' 'FC5' 'FC6' 'F7' 'F8'};
chans2use = [];
for chani = 1:length(channames)
    chans2use(chani) = find( strcmpi( channames{chani}, {EEG.chanlocs.labels} ));
end

powspec  = log10( mean( squeeze( mean( mean( mean( cleandatfr(:, :, chans2use, :, :) ,3,'omitnan' ) ,2,'omitnan' ) ,1,'omitnan' )  ) ,2,'omitnan') ) ;



% plot over time
powtopotheta30  = log10( squeeze( mean( mean( mean( mean( cleandatfr(:, :, :, thetaidx(1):thetaidx(2), 1:31)    ,5,'omitnan') ,4,'omitnan') ,2,'omitnan') ,1,'omitnan') ) );
powtopotheta150 = log10( squeeze( mean( mean( mean( mean( cleandatfr(:, :, :, thetaidx(1):thetaidx(2), 150:180) ,5,'omitnan') ,4,'omitnan') ,2,'omitnan') ,1,'omitnan') ) );
powtopotheta350 = log10( squeeze( mean( mean( mean( mean( cleandatfr(:, :, :, thetaidx(1):thetaidx(2), 320:350) ,5,'omitnan') ,4,'omitnan') ,2,'omitnan') ,1,'omitnan') ) );
powtopotheta500 = log10( squeeze( mean( mean( mean( mean( cleandatfr(:, :, :, thetaidx(1):thetaidx(2), 470:500) ,5,'omitnan') ,4,'omitnan') ,2,'omitnan') ,1,'omitnan') ) );

figure(5); colormap jet
subplot(1,4,1); topoplotIndie(powtopotheta30,  EEG.chanlocs(1:30),'numcontour',0,'electrodes','numbers','shading','interp');
title({'6.5 - 7.5 Hz,  epoch 1-31'});    c = colorbar; c.Label.String = 'Power \muV^2'; set(gca,'clim',[0.2 0.9]);
subplot(1,4,2); topoplotIndie(powtopotheta150, EEG.chanlocs(1:30),'numcontour',0,'electrodes','numbers','shading','interp');
title({'6.5 - 7.5  Hz,  epoch 150-180'}); c = colorbar; c.Label.String = 'Power \muV^2'; set(gca,'clim',[0.2 0.9]);
subplot(1,4,3); topoplotIndie(powtopotheta350, EEG.chanlocs(1:30),'numcontour',0,'electrodes','numbers','shading','interp');
title({'6.5 - 7.5 Hz,  epoch 320-350'}); c = colorbar; c.Label.String = 'Power \muV^2'; set(gca,'clim',[0.2 0.9]);
subplot(1,4,4); topoplotIndie(powtopotheta500, EEG.chanlocs(1:30),'numcontour',0,'electrodes','numbers','shading','interp');
title({'6.5 - 7.5 Hz,  epoch 470-500'}); c = colorbar; c.Label.String = 'Power \muV^2'; set(gca,'clim',[0.2 0.9]);

% Looks like overall (condition-general) theta amplitude is highest over
% Fz, and around these electrodes: {'AF3'}    {'FC1'}    {'Fz'}    {'Cz'}    {'AF4'}    {'FC2'}
% EEG chanloc indices: 2 4 29 30 27 25

conditions = readtable('Random_allocation_log_PITA_deblinded_matlab.csv');
subjects = table2array(conditions(:,1));
subjmatches = find( ismember(subjects, subj_list) > 0 );
conds = conditions(subjmatches,:);
s1_real = find( table2array(conds(:,2)) >  0 )';
s1_real= find( table2array(conds(:,2)) == 0 )';
s2_real = find( table2array(conds(:,3)) >  0 )';
s2_sham = find( table2array(conds(:,3)) == 0 )';

powtopotheta30_sham  = mean( [log10( squeeze( mean( mean( mean( cleandatfr(s1_sham, 1, :, thetaidx(1):thetaidx(2), 1:31)    ,5,'omitnan') ,4,'omitnan') ,1,'omitnan') ) ) ...
                              log10( squeeze( mean( mean( mean( cleandatfr(s2_sham, 2, :, thetaidx(1):thetaidx(2), 1:31)    ,5,'omitnan') ,4,'omitnan') ,1,'omitnan') ) )] ,2);
powtopotheta150_sham = mean( [log10( squeeze( mean( mean( mean( cleandatfr(s1_sham, 1, :, thetaidx(1):thetaidx(2), 150:180) ,5,'omitnan') ,4,'omitnan') ,1,'omitnan') ) ) ...
                              log10( squeeze( mean( mean( mean( cleandatfr(s2_sham, 2, :, thetaidx(1):thetaidx(2), 150:180) ,5,'omitnan') ,4,'omitnan') ,1,'omitnan') ) )] ,2);
powtopotheta350_sham = mean( [log10( squeeze( mean( mean( mean( cleandatfr(s1_sham, 1, :, thetaidx(1):thetaidx(2), 320:350) ,5,'omitnan') ,4,'omitnan') ,1,'omitnan') ) ) ...
                              log10( squeeze( mean( mean( mean( cleandatfr(s2_sham, 2, :, thetaidx(1):thetaidx(2), 320:350) ,5,'omitnan') ,4,'omitnan') ,1,'omitnan') ) )] ,2);
powtopotheta500_sham = mean( [log10( squeeze( mean( mean( mean( cleandatfr(s1_sham, 1, :, thetaidx(1):thetaidx(2), 470:500) ,5,'omitnan') ,4,'omitnan') ,1,'omitnan') ) ) ...
                              log10( squeeze( mean( mean( mean( cleandatfr(s2_sham, 2, :, thetaidx(1):thetaidx(2), 470:500) ,5,'omitnan') ,4,'omitnan') ,1,'omitnan') ) )] ,2);
powtopotheta30_real  = mean( [log10( squeeze( mean( mean( mean( cleandatfr(s1_real, 1, :, thetaidx(1):thetaidx(2), 1:31)    ,5,'omitnan') ,4,'omitnan') ,1,'omitnan') ) ) ...
                              log10( squeeze( mean( mean( mean( cleandatfr(s2_real, 2, :, thetaidx(1):thetaidx(2), 1:31)    ,5,'omitnan') ,4,'omitnan') ,1,'omitnan') ) )] ,2);
powtopotheta150_real = mean( [log10( squeeze( mean( mean( mean( cleandatfr(s1_real, 1, :, thetaidx(1):thetaidx(2), 150:180) ,5,'omitnan') ,4,'omitnan') ,1,'omitnan') ) ) ...
                              log10( squeeze( mean( mean( mean( cleandatfr(s2_real, 2, :, thetaidx(1):thetaidx(2), 150:180) ,5,'omitnan') ,4,'omitnan') ,1,'omitnan') ) )] ,2);
powtopotheta350_real = mean( [log10( squeeze( mean( mean( mean( cleandatfr(s1_real, 1, :, thetaidx(1):thetaidx(2), 320:350) ,5,'omitnan') ,4,'omitnan') ,1,'omitnan') ) ) ...
                              log10( squeeze( mean( mean( mean( cleandatfr(s2_real, 2, :, thetaidx(1):thetaidx(2), 320:350) ,5,'omitnan') ,4,'omitnan') ,1,'omitnan') ) )] ,2);
powtopotheta500_real = mean( [log10( squeeze( mean( mean( mean( cleandatfr(s1_real, 1, :, thetaidx(1):thetaidx(2), 470:500) ,5,'omitnan') ,4,'omitnan') ,1,'omitnan') ) ) ...
                              log10( squeeze( mean( mean( mean( cleandatfr(s2_real, 2, :, thetaidx(1):thetaidx(2), 470:500) ,5,'omitnan') ,4,'omitnan') ,1,'omitnan') ) )] ,2);


figure(7); colormap jet
subplot(2,4,1); topoplotIndie(powtopotheta30_sham,  EEG.chanlocs(1:30),'numcontour',0,'electrodes','numbers','shading','interp');
title({'SHAM, 4.5-5.5 Hz, ep:1-31'});    c = colorbar; c.Label.String = 'Power \muV^2'; set(gca,'clim',[0.2 0.9]);
subplot(2,4,2); topoplotIndie(powtopotheta150_sham, EEG.chanlocs(1:30),'numcontour',0,'electrodes','numbers','shading','interp');
title({'SHAM, 4.5-5.5 Hz, ep:150-180'}); c = colorbar; c.Label.String = 'Power \muV^2'; set(gca,'clim',[0.2 0.9]);
subplot(2,4,3); topoplotIndie(powtopotheta350_sham, EEG.chanlocs(1:30),'numcontour',0,'electrodes','numbers','shading','interp');
title({'SHAM, 4.5-5.5 Hz, ep:320-350'}); c = colorbar; c.Label.String = 'Power \muV^2'; set(gca,'clim',[0.2 0.9]);
subplot(2,4,4); topoplotIndie(powtopotheta500_sham, EEG.chanlocs(1:30),'numcontour',0,'electrodes','numbers','shading','interp');
title({'SHAM, 4.5-5.5 Hz, ep:470-500'}); c = colorbar; c.Label.String = 'Power \muV^2'; set(gca,'clim',[0.2 0.9]);
subplot(2,4,5); topoplotIndie(powtopotheta30_real,  EEG.chanlocs(1:30),'numcontour',0,'electrodes','numbers','shading','interp');
title({'REAL, 4.5-5.5 Hz, ep:1-31'});    c = colorbar; c.Label.String = 'Power \muV^2'; set(gca,'clim',[0.2 0.9]);
subplot(2,4,6); topoplotIndie(powtopotheta150_real, EEG.chanlocs(1:30),'numcontour',0,'electrodes','numbers','shading','interp');
title({'REAL, 4.5-5.5 Hz, ep:150-180'}); c = colorbar; c.Label.String = 'Power \muV^2'; set(gca,'clim',[0.2 0.9]);
subplot(2,4,7); topoplotIndie(powtopotheta350_real, EEG.chanlocs(1:30),'numcontour',0,'electrodes','numbers','shading','interp');
title({'REAL, 4.5-5.5 Hz, ep:320-350'}); c = colorbar; c.Label.String = 'Power \muV^2'; set(gca,'clim',[0.2 0.9]);
subplot(2,4,8); topoplotIndie(powtopotheta500_real, EEG.chanlocs(1:30),'numcontour',0,'electrodes','numbers','shading','interp');
title({'REAL, 4.5-5.5 Hz, ep:470-500'}); c = colorbar; c.Label.String = 'Power \muV^2'; set(gca,'clim',[0.2 0.9]);

mean30_sham  = mean( powtopotheta30_sham(  [2 4 29 30 27 25] ) );
std30_sham   = std(  powtopotheta30_sham(  [2 4 29 30 27 25] ) );
mean150_sham = mean( powtopotheta150_sham( [2 4 29 30 27 25] ) );
std150_sham  = std(  powtopotheta150_sham( [2 4 29 30 27 25] ) );
mean350_sham = mean( powtopotheta350_sham( [2 4 29 30 27 25] ) );
std350_sham  = std(  powtopotheta350_sham( [2 4 29 30 27 25] ) );
mean500_sham = mean( powtopotheta500_sham( [2 4 29 30 27 25] ) );
std500_sham  = std(  powtopotheta500_sham( [2 4 29 30 27 25] ) );
mean30_real  = mean( powtopotheta30_real(  [2 4 29 30 27 25] ) );
std30_real   = std(  powtopotheta30_real(  [2 4 29 30 27 25] ) );
mean150_real = mean( powtopotheta150_real( [2 4 29 30 27 25] ) );
std150_real  = std(  powtopotheta150_real( [2 4 29 30 27 25] ) );
mean350_real = mean( powtopotheta350_real( [2 4 29 30 27 25] ) );
std350_real  = std(  powtopotheta350_real( [2 4 29 30 27 25] ) );
mean500_real = mean( powtopotheta500_real( [2 4 29 30 27 25] ) );
std500_real  = std(  powtopotheta500_real( [2 4 29 30 27 25] ) );

numbers_theta = [mean30_sham  mean30_real std30_sham std30_real; ...
    mean150_sham  mean150_real std150_sham std150_real; ...
    mean350_sham  mean350_real std350_sham std350_real; ...
    mean500_sham  mean500_real std500_sham std500_real];

timeaxis2 = [30 150 350 500];

errlow  = numbers_theta(:,1:2) - numbers_theta(:,3:4);
errhigh = numbers_theta(:,1:2) + numbers_theta(:,3:4);

figure(10)
bar( numbers_theta(:,1:2))

hold on

er = errorbar( numbers_theta(:,1:2),errlow,errhigh);    
%er.Color = [0 0 0];                            
%er.LineStyle = 'none';  

hold off


figure(10)

% data 
model_series = [mean30_sham  mean30_real ; ...
    mean150_sham  mean150_real ; ...
    mean350_sham  mean350_real ; ...
    mean500_sham  mean500_real ];
model_error = [std30_sham std30_real; ...
    std150_sham std150_real; ...
    std350_sham std350_real; ...
    std500_sham std500_real];
b = bar(model_series, 'grouped');
hold on
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(model_series);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',model_series,model_error,'k','linestyle','none');
hold off

[h, p] = ttest(x,y)

figure(10)
plot([30 150 350 500], model_series(:,1), 'k', 'LineWidth', 4 ) ; set(gca,'ylim',[0.5 0.8])
hold on
plot([30 150 350 500], model_series(:,2), 'r', 'LineWidth', 4 ) ;
% ylabel({'Power muV^2')

% frqbands IFCN 2020 guidelines (TF = transition freq between theta and alpha, IAF = individual alpha freq peak):
% delta: 0.1-<4 Hz |or| from TF−4Hz to TF−2Hz
% theta: 4-<8   Hz |or| from TF−2Hz to TF  
% alpha: 8-13   Hz |or| from TF to IAF+2Hz 
% beta:  14-30  Hz 
% gamma: >30-80 Hz

% frequency boundaries in Hz & convert to indices ()
thetawin  = dsearchn(hz,[4  8]');
alphawin  = dsearchn(hz,[8 13]');

% Individual alpha peak frequency
alphapows = log10( mean( powspec(chani, alphawin(1):alphawin(2), : ) ,3 ) ); %mean powerspectrum for alpha frequencies in this channel, 1/F remove by logtransform
figure(2); plot( hz(alphawin(1):alphawin(2)), alphapows )
IAFidx    = find( alphapows == max(alphapows) );
IAF       = alphawin(1) + IAFidx - 1;
