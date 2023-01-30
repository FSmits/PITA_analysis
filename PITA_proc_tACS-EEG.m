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




% %% Verify stimulator was on during tACS procedure in each individual dataset
% 
% % enter subject identification numbers
% full_subj_list =   [669	557 363	638	602	989	383	502	733	442	575	913 710	262 ...
%     752 227	565	626 334	362	600	121	319 923	915	298	202	692	275	...
%     508 291	803	755	681	876	134	559	396	818	601	297	524	883	193	642];
% 
% for subj_i = 1:length(full_subj_list)
%     for sess_i = 1:length(sessions)
% 
%         fprintf('\n****\nStart processing subject %i session %i\n****\n\n', full_subj_list(subj_i), sessions(sess_i));
% 
%         fileName = ['TACSEEG-' num2str(full_subj_list(subj_i)) '-' num2str(sess_i)];
% 
%         % Load EEG set
%         EEG      = pop_loadset('filename', [fileName, '_RawEEG.set'], 'filepath', Path2EEGsets);
% 
%         % Plot data and check if tACS artifact is visible. If not - check if stimulation was carried out correctly.
%         pop_eegplot( EEG, 1, 1, 1);
% 
%         % Check if tACS has been carried out by checking the typical tACS artifact
%         m = 0 ;
%         while m == 0
%             m = input('Do you see tACS artifacts? Y/N:','s');
%             if m == 'N'
%                 writecell( { ['tACS-EEG-' num2str(full_subj_list(subj_i)) '-' num2str(sess_i)] } , ...
%                     ['/Users/fsmits2/Downloads/1 EEG data tacs-eeg/3 EEG sets processed/' ...
%                     ['noStim-' num2str(full_subj_list(subj_i)) '-' num2str(sess_i) '.txt']]);
%                 continue
%             elseif m == 'Y'
%                 counter = counter + 1;
%             end
%         end
%     end
% end

%% Start by clearing workspace

clear
close all


%% Initialize

% open EEGlab
cd('/Users/fsmits2/Downloads/eeglab2022.1')
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

cd('/Users/fsmits2/Documents/PITA_analysis'); % return to PITA analysis folder

% set paths to data
Path2EEGsets = '/Users/fsmits2/Downloads/1 EEG data tacs-eeg/3 EEG sets processed';

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

bdchns        = table2cell(  readtable( [Path2EEGsets '/Overview_badchannels_20-Jan-2023.txt'] ,'Format','auto') );




subj_i = 1
sess_i = 2

        fprintf('\n****\nLoad tACS-EEG preprocessed data from: subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));
        fileName = [file_type{3} num2str(subj_list(subj_i)) '-' num2str(sess_i) '_CleanEEG.set'];

        % Load EEG set
        EEG      = pop_loadset('filename', fileName, 'filepath', Path2EEGsets);



         % Compute the laplacian 
    EEG = pop_chanedit(EEG, 'lookup','H:\Myrte\eeglab2022.0\plugins\dipfit\standard_BESA\standard-10-5-cap385.elp');
    
    cd('H:\Myrte\Scripts, Functies & Output'); 
    EEG.data = laplacian_perrinX(EEG.data,[EEG.chanlocs.X],[EEG.chanlocs.Y],[EEG.chanlocs.Z]);

        % ------ Do the FFT -------



       % ----- Select only from channels with clean data, i.e. remove noisy channels from analysis outcome:

        badchannels = bdchns{subj_i, sess_i+1};
        chanarray   = 1:length(EEG.chanlocs); 
        if sum( strcmpi( badchannels, '0') ) < 1
        badchans = regexp(badchannels, ',', 'split')
        noisychannels = [];
            for bchni = 1:length(badchans)
                noisychannels(bchni) = find( strcmpi( badchans{bchni}, {EEG.chanlocs.labels} ));
            end
        chanarray(noisychannels) = []; % remove the noisy channels from chanarray
        end


