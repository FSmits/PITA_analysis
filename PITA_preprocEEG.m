%% EEG pre-processing - PITA - all data

clear
close all

%% open EEGlab 

cd('/Users/fsmits2/Downloads/eeglab2021.0')
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

% return to PITA analysis folder
cd('/Users/fsmits2/Documents/PITA_analysis');


%% set paths to data

Path2EEGbdf  = '/Users/fsmits2/Downloads/';
Path2EEGsets = '/Users/fsmits2/Downloads';


%% enter subject identification numbers

subj = [669	557	363	638	602	989	383	502	733	442	575	913	710	262 ...
    752	227	565	626	334	362	600	121	319	923	915	298	202	692	275	...
    508 291	803	755	681	876	134	559	396	818	601	297	524	883	193	642];

%% enter filenames

rec1 = 'restingstate_EEG_pretACS';

