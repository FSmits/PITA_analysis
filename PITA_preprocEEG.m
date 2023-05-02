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

% set paths to data
Path2EEGbdf  = '/Users/fsmits2/Downloads/1 EEG data tacs-eeg/';
Path2EEGsets = '/Users/fsmits2/Downloads/1 EEG data tacs-eeg/post-tacs r-s processed';
path2save    = '/Users/fsmits2/Documents/PITA_analysis';
path2scripts = '/Users/fsmits2/MATLAB/Projects/PITA_analysissandbox';

cd(path2scripts)

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


%% read EEG file

% Original trigger codes
StartRec        = 254;
StopRec         = 255;

% resting-state EEG
EyesOpenOnset    = 1;
EyesOpenOffset   = 2;
EyesClosedOnset  = 3;
EyesClosedOffset = 4;

% tACS-EEG
start_tACS      = 5; %(2:20) identical to stop_EEG minus(1:19)
stop_tACS       = 6; %(1:20) identical to start_EEG (1:20)
start_EEG       = 7;
stop_EEG        = 8;
post_tACS_EEG   = 9;

% %%% For tACS-EEG data:
% % Initiate time/period-related triggers
% secs      = 30; % 30 seconden data na elke tACS stimulatie
% stims     = 20; % 20 tACS stimulaties in totaal
% trig_base = repmat(0.01:0.01:0.30,[stims,1]);
% stim_mat  = repmat(1:stims,[secs,1])';
% trigs     = trig_base + stim_mat; % Trigger names are: #stimulation as integer, #second of data following that stimulation as decimal

% Which task (file type) you want to analyze?
fileno = 3;

% Loop over files
for subj_i = 31:length(subj_list)
    for sess_i = 1:length(sessions)

        fprintf('\n****\nStart processing subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));
        fileName = fullfile(Path2EEGbdf, [file_type{fileno} num2str(subj_list(subj_i)) '-' num2str(sess_i) '.bdf']);

        % -- Load raw bdf data file via EEGlab
        EEG = pop_biosig( fileName );

        % -- Enter data to the EEG structure
        EEG.filename = fileName;
        EEG.setname  = fileName;
        EEG.subject  = subj_list(subj_i);
        EEG.session  = sess_i;

        % -- Remove non-recorded channels: F3 F4 (tACS electrode locations) and EXG7 EXG8
        EEG = pop_select(EEG, 'nochannel', {'F3', 'F4', 'EXG7', 'EXG8'});

        % -- Re-code events [Why? BioSemi/Computer settings resulted in changes in the recorded trigger codes relative to the originally programmed triggers codes. These changes are unfortunately not exactly the same across subjects.]
        % remove added trigger text like 'condition' and 'artifact'
        for ev_i = 1:length({EEG.event.type})
            EEG.event(ev_i).type = strrep( EEG.event(ev_i).type, 'condition ', '' );
            EEG.event(ev_i).type = strrep( EEG.event(ev_i).type, 'artifact', '' );
        end

        % remove trigger 256
        trig_256 = find(strcmpi( {EEG.event.type}, '256' ));
        EEG = pop_editeventvals(EEG,'delete', trig_256);


%         %%% For pre-tACS resting-state data:
%         % Delete double events
%         if subj_list(subj_i)==710 && sess_i==1
%             EEG = pop_editeventvals(EEG,'delete', [3 11]);
%         elseif subj_list(subj_i)==752 && sess_i==1
%             EEG = pop_editeventvals(EEG,'delete', [6 8 13]);
%         elseif subj_list(subj_i)==752 && sess_i==2
%             EEG = pop_editeventvals(EEG,'delete', [1 2 3]);
%         elseif subj_list(subj_i)==600 && sess_i==2
%             EEG = pop_editeventvals(EEG,'delete', [5 9]);
%         elseif subj_list(subj_i)==923 && sess_i==2
%             EEG = pop_editeventvals(EEG,'delete', [4]);
%         elseif subj_list(subj_i)==508 && sess_i==2
%             EEG = pop_editeventvals(EEG,'delete', [3 6]);
%             EEG = pop_editeventvals(EEG,'insert',{ 4 ,[],[],[],[]},'changefield',{ 4 ,'type','3' }, 'changefield',{ 4 ,'edftype','3' }, 'changefield',{ 4 ,'latency', EEG.event(3).latency/EEG.srate+1 }); % latency of events (EEG.event.latency) is defined in data points, not in time. But when you want to change it, you need to define in seconds.
%         elseif subj_list(subj_i)==681 && sess_i==2
%             EEG = pop_editeventvals(EEG,'delete', [1 2 8 13]);
%         elseif subj_list(subj_i)==876 && sess_i==2
%             EEG = pop_editeventvals(EEG,'delete', [8 10]);
%         elseif subj_list(subj_i)==134 && sess_i==2
%             EEG = pop_editeventvals(EEG,'delete', [3]);
%         elseif subj_list(subj_i)==502 && sess_i==2
%             EEG = pop_editeventvals(EEG,'delete', [9]);
%             EEG = pop_editeventvals(EEG,'insert',{ 2 ,[],[],[],[]},'changefield',{ 2 ,'type','2' }, 'changefield',{ 2 ,'edftype','2' }, 'changefield',{ 2 ,'latency', EEG.event(3).latency/EEG.srate-1 }); % latency of events (EEG.event.latency) is defined in data points, not in time. But when you want to change it, you need to define in seconds.
%         end

        %%% For post-tACS resting-state data:
        if fileno==3 %when post-tACS resting-state is inside tACS-EEG bdf datatset
            rs_begin = [find(strcmpi( {EEG.event.type}, '254' )) find(strcmpi( {EEG.event.type}, '510' ))];
            if subj_list(subj_i)==681 && sess_i==2
                rs_begin = rs_begin(1);
            end
            if length(rs_begin)>1 || isempty(rs_begin)
                fprintf('****\nNone or too many resting-state start triggers in subject %i session %i\n', subj_list(subj_i), sessions(sess_i));
                return
            end
            EEG = eeg_eegrej(EEG, [0,  EEG.event(rs_begin).latency]); % remove data during tACS
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG); % store changes
        end

        % Delete double events
        if subj_list(subj_i)==600 && sess_i==1 && fileno==4
            EEG = pop_editeventvals(EEG,'delete', [9]);
        end

        % find or create trigger for resting-state onset  eyes-open & eyes-closed  (original codes = '1' / '3'. Other possibile codes = ':EOG (blinks, fast, large amplitude)' / ':EMG/Muscle'.)
        open_begin  = sort([find(strcmpi( {EEG.event.type}, '1' )), ...
            find(strcmpi( {EEG.event.type}, ':EOG (blinks, fast, large amplitude)' )) ]);
        close_begin = sort([find(strcmpi( {EEG.event.type}, '3' )), ...
            find(strcmpi( {EEG.event.type}, ':EMG/Muscle' )) ]);     

        % find trigger for eyes-open & eyes-closed resting-state offset  (original code = '2' / '4'. Other possibile codes = ':ECG' / ':Movement'.)
        open_end  = sort([find(strcmpi( {EEG.event.type}, '2' )), ...
            find(strcmpi( {EEG.event.type}, ':ECG' )) ]);
        close_end = sort([find(strcmpi( {EEG.event.type}, '4' )), ...
            find(strcmpi( {EEG.event.type}, ':Movement' )) ]);

        % re-code to resting-state onsets and offsets
        for trig_i = 1:length(open_begin)
            EEG = pop_editeventvals(EEG,'changefield', {open_begin(trig_i)  'type' 'EyesOpenOnset'});
        end
        for trig_i = 1:length(open_end)
            EEG = pop_editeventvals(EEG,'changefield', {open_end(trig_i)    'type' 'EyesOpenOffset'});
        end
        for trig_i = 1:length(close_begin)
            EEG = pop_editeventvals(EEG,'changefield', {close_begin(trig_i) 'type' 'EyesClosedOnset'});
        end
        for trig_i = 1:length(close_end)
            EEG = pop_editeventvals(EEG,'changefield', {close_end(trig_i)   'type' 'EyesClosedOffset'});
        end
  

% % %         %%% For tACS-EEG data:
% % %         % find or create trigger for tACS onset   (original code = '5'. Other possibile codes = '8' , ':Failing electrode' , ':breathing'.)
% % %         if any( subj_list(subj_i) == [638 275 262] )  &&  sess_i==1
% % %             EEG = pop_editeventvals(EEG,'insert',{1,[],[],[],[]},'changefield',{1,'type','5'},'changefield',{1,'edftype','5'},'changefield',{1,'latency',0.01});
% % %         end
% % %         tACS_begin = sort([find(strcmpi( {EEG.event.type}, '5' )), ...
% % %             find(strcmpi( {EEG.event.type}, '8' )),...
% % %             find(strcmpi( {EEG.event.type}, ':Failing electrode' )),...
% % %             find(strcmpi( {EEG.event.type}, ':breathing' )) ]);
% % % 
% % %         % find trigger for tACS offset   (original code = '6'. Other possibile codes = '7' , ':sweat' , ':50/60 Hz mains interference'.)
% % %         tACS_end = sort([find(strcmpi( {EEG.event.type}, '6' )), ...
% % %             find(strcmpi( {EEG.event.type}, '7' )), ...
% % %             find(strcmpi( {EEG.event.type}, ':50/60 Hz mains interference' )), ...
% % %             find(strcmpi( {EEG.event.type}, ':sweat' )) ]);
% % % 
% % %         % when too many trigger codes are found, select codes that show >90-second trigger-to-trigger delay (tACS_begin and tACS_end triggers are separated by ~92 seconds: 60 seconds stimulation + 30 second EEG recording period + 2 seconds ramp-up/ramp-down stimulation)
% % %         if length(tACS_begin) > 21
% % %             evlat90_begin = [];
% % %             evlat = [];
% % %             for eventi = 1:length(tACS_begin)
% % %                 evlat(eventi) = EEG.event(tACS_begin(eventi)).latency / EEG.srate; %save event latencies in array
% % %             end
% % %             evlat90_begin  = find( diff(evlat) > 90) + 1; %find which events have a latency difference of >90 seconds
% % %             tACS_begin     = [tACS_begin(1) tACS_begin(evlat90_begin)];
% % %         end
% % %         if length(tACS_end) > 21
% % %             evlat90_end = [];
% % %             evlat = [];
% % %             for eventi = 1:length(tACS_end)
% % %                 evlat(eventi) = EEG.event(tACS_end(eventi)).latency / EEG.srate; %save event latencies in array
% % %             end
% % %             evlat90_end  = find( diff(evlat) > 90) + 1; %find which events have a latency difference of >90 seconds
% % %             tACS_end     = [tACS_end(1) tACS_end(evlat90_end)];
% % %         end
% % % 
% % %         % re-code to original code for tACS onsets and offsets
% % %         for trig_i = 1:length(tACS_end)
% % %             EEG = pop_editeventvals(EEG,'changefield', {tACS_begin(trig_i) 'type' 'tACS_start'});
% % %             EEG = pop_editeventvals(EEG,'changefield', {tACS_end(trig_i)   'type' 'tACS_stop'});
% % %         end
% % % 
% % %         % insert time-related indices as event codes
% % %         bndrs     = find(strcmpi( {EEG.event.type}, 'tACS_stop' ));
% % %         bndrs_lat = [EEG.event(bndrs).latency];
% % %         % insert indexing trigger code every full second + little delay of 5 ms so that the indexing triggers do not overlap with 1-second epoch cuts
% % %         for i_EEGs = 1:length(bndrs)
% % %             for i_secs = 1:secs
% % %                 EEG = pop_editeventvals(EEG, 'insert',{bndrs(i_EEGs)+i_secs,[],[],[],[]},...
% % %                     'changefield',{bndrs(i_EEGs)+i_secs ,'type',    trigs(i_EEGs,i_secs) },...
% % %                     'changefield',{bndrs(i_EEGs)+i_secs ,'edftype', trigs(i_EEGs,i_secs) },...
% % %                     'changefield',{bndrs(i_EEGs)+i_secs ,'latency', bndrs_lat(i_EEGs)/EEG.srate + i_secs}); % for latency event, latencies are in millisecond compared to the time locking event, not in data samples.
% % %             end
% % %         end
% % % 
% % %         % -- Cut out data during stimulation (tACS)
% % %         tACS_begin = sort(find(strcmpi( {EEG.event.type}, 'tACS_start' )));
% % %         tACS_end   = sort(find(strcmpi( {EEG.event.type}, 'tACS_stop' )));
% % %         for ni = flip(1:length(tACS_begin))
% % %             EEG = eeg_eegrej(EEG, [EEG.event(tACS_begin(ni)).latency,  EEG.event(tACS_end(ni)).latency]); % remove data during tACS
% % %             [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG); % store changes
% % %         end
% % % 
% % %         % -- Cut out data before first tACS-EEG & after last tACS-EEG period
% % %         EEG_begin = find(strcmpi( {EEG.event.type}, 'boundary' ));
% % %         EEG_end   = sort([find(strcmpi( {EEG.event.type}, '8' )), find(strcmpi( {EEG.event.type}, ':breathing' ))]);
% % %         EEG_end   = EEG_end(1,end); % in case multiple trigger codes '8' are available, pick last one
% % %         EEG       = eeg_eegrej(EEG, [EEG.event(EEG_end).latency + 0.05,   EEG.pnts(end)]);
% % %         EEG       = eeg_eegrej(EEG, [0,   EEG.event(EEG_begin(1)).latency - 0.05]);

        %         pop_eegplot( EEG, 1, 1, 1); % Inspect data

        % -- Save
        fprintf('\n****\nSave processed subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));
        SaveName = [file_type{fileno} num2str(subj_list(subj_i)) '-' num2str(sess_i) '_RawEEG.set'];
        EEG = pop_saveset( EEG, 'filename',SaveName,'filepath', Path2EEGsets );

        clear EEG
        ALLEEG(1:end) = [];
    end
end



%% Pre-processing steps: re-reference, downsample, filter, create bipolar EOG channels

fileno = 4;

for subj_i = 11:length(subj_list)
    for sess_i = 1:length(sessions)

        fprintf('\n****\nStart pre-processing subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));
        fileName = [file_type{fileno} num2str(subj_list(subj_i)) '-' num2str(sess_i) '_RawEEG.set'];

        % Load EEG set
        EEG      = pop_loadset('filename', fileName, 'filepath', Path2EEGsets);

% % %         %%% For tACS-EEG
% % %         % Check if all EEG-recording episodes are present
% % %         m0 = -1 ;
% % %         fprintf([ num2str(subj_list(subj_i)) ' session ' num2str(sessions(sess_i)) ' -- Number of events: ' num2str(length(EEG.event)) ' \n '])
% % %         while m0 == -1
% % %             if length(EEG.event) < 600
% % %                 m0 = input('Continue? Y/N: ','s');
% % %                 if m0 == 'Y'
% % %                     continue
% % %                 else 
% % %                     m0 = 0;
% % %                     return
% % %                 end
% % %             else
% % %                 m0 = 0;
% % %                 continue
% % %             end
% % %         end

        % Re-reference to avg mastoids
        mastoid1 = find(strcmpi( {EEG.chanlocs.labels}, 'EXG5' ));
        mastoid2 = find(strcmpi( {EEG.chanlocs.labels}, 'EXG6' ));
        EEG      = pop_reref( EEG, [mastoid1 mastoid2]); %re-references to the average of 2 channels

        % Downsample & Filter
        EEG      = pop_resample( EEG, 256); % Downsample the data from 2048 to 256 Hz
        EEG      = pop_basicfilter( EEG, 1:32 , 'Cutoff',  0.5, 'Design', 'butter', 'Filter', 'highpass', 'Order',  4 ); % Format: pop_basicfilter( EEG, chanArray, parameters )
        EEG      = pop_basicfilter( EEG, 1:32 , 'Cutoff',   35, 'Design', 'butter', 'Filter',  'lowpass', 'Order',  4 ); % IIR Butterworth filters highpass 0.5 Hz, lowpass 35 Hz, filter order 4 (-24 dB rolloff).

        % Create eye channel bipolar signals
        EOG_ch   = find(strcmpi({EEG.chanlocs.labels},'EXG1')):find(strcmpi({EEG.chanlocs.labels},'EXG4')); % Find indices EOG electrodes (EXG1, EXG2, EXG3, EXG4) & EEG scalp electrodes
        EEG_ch   = 1:EOG_ch(1)-1;
        EEG      = pop_reref(EEG, EOG_ch(2),'exclude',[EEG_ch, EOG_ch(3:4)] ); % For vertical eye movements: re-reference channel below left eye (EXG1) to channel above left eye (EXG2):
        EOG_ch   = find(strcmpi( {EEG.chanlocs.labels},'EXG1')):find(strcmpi( {EEG.chanlocs.labels},'EXG4')); % Find indices EOG electrodes (EXG1, EXG2, EXG3, EXG4) & EEG scalp electrodes
        EEG      = pop_reref(EEG, EOG_ch(2),'exclude',[EEG_ch, EOG_ch(1)] ); % For horizontal eye movements: re-reference channel next to left eye. (EXG3) to channel next to right eye (EXG4):
        EEG.chanlocs( EOG_ch(1) ).labels = 'VEOG'; % EXG1 is now the bipolar VEOG channel. Change channel name.
        EEG.chanlocs( EOG_ch(2) ).labels = 'HEOG'; % EXG3 is now the bipolar HEOG channel. Change channel name.

%        % CleanLine (Notch filter via ICA) to remove 50 Hz power line noise
%        EEG      = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:30] ,'computepower',1,'linefreqs',50,'newversion',0,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',0,'sigtype','Channels','taperbandwidth',2,'tau',100,'verb',1,'winsize',4,'winstep',1);

        %         pop_eegplot( EEG, 1, 1, 1); % Inspect data

        % Save
        fprintf('\n****\nSave pre-processed subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));
        SaveName = [file_type{fileno} num2str(subj_list(subj_i)) '-' num2str(sess_i) '_PreprocEEG.set'];
        EEG = pop_saveset( EEG, 'filename',SaveName,'filepath', Path2EEGsets );

        clear EEG
        ALLEEG(1:end) = [];
    end
end



%% Interpolate channels and detect gel bridges

% % Loop over files
% for subj_i = 1:length(subj_list)
%     for sess_i = 1:length(sessions)
%         k0 = [];
%         k0 = input('Gel bridge? [1 = yes / 0 = no] ');
%         if k0 == 1
%             gelbridges(subj_i,sess_i+1) = 1;
%         else
%             gelbridges(subj_i,sess_i+1) = 0;
%         end
% 
%         writematrix(gelbridges, [Path2EEGsets '/Overview_gelbridges_'           char(datetime('today')) '.txt'], 'Delimiter',';'); %char(datetime('yesterday')) '.txt'], 'Delimiter',';'); %
%     end
% end


%% ICA

fileno = 4;

% Loop over files
for subj_i = 35:length(subj_list)
    for sess_i = 1:length(sessions)

        fprintf('\n****\nLoad subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));
        fileName = [file_type{fileno} num2str(subj_list(subj_i)) '-' num2str(sess_i) '_PreprocEEG.set'];

        % Load EEG set
        EEG      = pop_loadset('filename', fileName, 'filepath', Path2EEGsets);

        % Find Eyeblink components with ICA (runica):
        EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');

        % Save
        fprintf('\n****\nSave ICA data subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));
        SaveName = [file_type{fileno} num2str(subj_list(subj_i)) '-' num2str(sess_i) '_Preproc_ICAfull.set'];
        EEG      = pop_saveset( EEG, 'filename',SaveName,'filepath', Path2EEGsets );
          
        clear EEG
        ALLEEG(1:end) = [];

    end
end


%% [resting-state only] Epoch the data

fileno = 4;

for subj_i = 34:length(subj_list)
    for sess_i = 1:length(sessions)

        fprintf('\n****\nLoad subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));
        fileName = [file_type{fileno} num2str(subj_list(subj_i)) '-' num2str(sess_i) '_Preproc_ICAfull.set'];

        % Load EEG set
        EEG      = pop_loadset('filename', fileName , 'filepath', Path2EEGsets);

        % Check (& fix) if enough and the right triggers are present
% % %        % exceptions for pre-tACS resting-state eeg:
% % %         if subj_list(subj_i)==502 && sess_i == 2
% % %             EEG = pop_editeventvals(EEG,'delete', length(EEG.event)-1);
% % %             EEG = pop_editeventvals(EEG,'insert',{ length(EEG.event)-1 ,[],[],[],[]}, 'changefield',{ length(EEG.event)-1 ,'type','EyesClosedOffset' },  'changefield',{ length(EEG.event)-1 ,'edftype','EyesClosedOffset' }, 'changefield',{ length(EEG.event)-1 ,'latency', EEG.event(end).latency/EEG.srate-1 }); % latency of events (EEG.event.latency) is defined in data points, not in time. But when you want to change it, you need to define in seconds.
% % %         elseif (subj_list(subj_i)==227 || subj_list(subj_i)==565) && sess_i == 2
% % %             EEG = pop_editeventvals(EEG,'insert',{ 1 ,[],[],[],[]}, 'changefield',{ 1 ,'type','EyesOpenOnset' },  'changefield',{ 1 ,'edftype','EyesOpenOnset' }, 'changefield',{ 1 ,'latency', (EEG.times(1)/1000)+1 });
% % %         elseif subj_list(subj_i)==915 && sess_i == 2
% % %             EEG = pop_editeventvals(EEG,'insert',{ length(EEG.event) ,[],[],[],[]}, 'changefield',{ length(EEG.event) ,'type','EyesClosedOffset' },  'changefield',{ length(EEG.event) ,'edftype','EyesClosedOffset' }, 'changefield',{ length(EEG.event) ,'latency',  (EEG.times(end)/1000)-1 });
% % %         elseif subj_list(subj_i)==202 && sess_i == 2
% % %             trig_20s = [find(strcmpi( {EEG.event.type}, '20' )) find(strcmpi( {EEG.event.type}, '21' )) find(strcmpi( {EEG.event.type}, '22' )) find(strcmpi( {EEG.event.type}, '23' )) find(strcmpi( {EEG.event.type}, '20' ))];
% % %             EEG = pop_editeventvals(EEG,'delete', trig_20s);
% % %         end
        % exceptions for post-tACS eeg sets:
        if subj_list(subj_i)==669 && sess_i == 2
            EEG = pop_editeventvals(EEG,'delete', 2);
        elseif subj_list(subj_i)==638 && sess_i == 2
            EEG = pop_editeventvals(EEG,'delete', [2 9]);
        elseif subj_list(subj_i)==989 && sess_i == 1
            EEG = pop_editeventvals(EEG,'delete', 6);
        elseif subj_list(subj_i)==227 && sess_i == 2
            EEG = pop_editeventvals(EEG,'delete', [2 4 8]);
        elseif subj_list(subj_i)==362 && sess_i == 2
            EEG = pop_editeventvals(EEG,'insert',{5,[],[],[],[],[]},'changefield',{5,'latency',122.94},'changefield',{5,'type','EyesClosedOffset'},'changefield',{5,'edftype',271});
            EEG = pop_editeventvals(EEG,'delete', [8]);
        elseif subj_list(subj_i)==298 && sess_i == 2
            EEG = pop_editeventvals(EEG,'delete', [7]);
        elseif subj_list(subj_i)==681 && sess_i == 2
            EEG = pop_editeventvals(EEG,'delete', [2]);
        elseif subj_list(subj_i)==559 && sess_i == 1
            EEG = pop_editeventvals(EEG,'insert',{7,[],[],[],[],[]},'changefield',{7,'latency',184.02},'changefield',{7,'type','EyesOpenOffset'},'changefield',{7,'edftype',322});

        end

        open_begin  = sort(find(strcmpi( {EEG.event.type}, 'EyesOpenOnset' ))); %Find the open/closed eyes onset/offset triggers
        close_end   = sort(find(strcmpi( {EEG.event.type}, 'EyesClosedOffset' )));
        if length(open_begin) ~= 2 || length(close_end) ~= 2
            if subj_list(subj_i)~=383 && sess_i~=1
                fprintf('\n****\nNot enough trigger codes in subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));
                return
            end
        end

        % To enter trigger after each second in continuous data, first find timings in data
        ntrg   = [open_begin  close_end];
        for trgi = 1:length(ntrg)
            open_begin  = sort(find(strcmpi( {EEG.event.type}, 'EyesOpenOnset' ))); %Find the open/closed eyes onset/offset triggers
            open_end    = sort(find(strcmpi( {EEG.event.type}, 'EyesOpenOffset' )));
            close_begin = sort(find(strcmpi( {EEG.event.type}, 'EyesClosedOnset' )));
            close_end   = sort(find(strcmpi( {EEG.event.type}, 'EyesClosedOffset' )));
            trgT0   = [open_begin  close_begin];
            trgTend = [open_end    close_end  ];
            T0   = EEG.event(trgT0(  trgi)).latency;
            if subj_list(subj_i)==383 && sess_i==1
                if trgi==2
                    Tend = (EEG.times(end)/1000)*256;
                elseif trgi==3
                    Tend = EEG.event(trgTend(2)).latency;
                else
                    Tend = EEG.event(trgTend(trgi)).latency;
                end
            elseif subj_list(subj_i)==710 && sess_i==1 && trgi==4
                Tend = (EEG.times(end)/1000)*EEG.srate;
            else
                Tend = EEG.event(trgTend(trgi)).latency;
            end
            % Create triggers each second (1-sec segments) minus the final second (because no full second left):
            trggLats = T0:EEG.srate:Tend-1*EEG.srate;
            for segi = 1:length(trggLats)
                if ismember(trgT0(trgi),open_begin)
                    trgcode = '11';
                elseif ismember(trgT0(trgi),close_begin)
                    trgcode = '22';
                else
                    fprintf('\n****\nNot able to match trigger codes in subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));
                    return
                end
                EEG = pop_editeventvals(EEG,'insert',{ segi ,[],[],[],[]},...
                    'changefield',{ segi ,'type',    trgcode },...
                    'changefield',{ segi ,'edftype', trgcode },...
                    'changefield',{ segi ,'latency', trggLats(segi)/256 }); % latency of events (EEG.event.latency) is defined in data points, not in time. But when you want to change it, you need to define in seconds.
            end
        end

        % Divide into 1-sec epochs
        EEG = pop_epoch( EEG, {'11' '22'},  [0  1], 'epochinfo', 'yes');
        % Remove baseline mean
        EEG = pop_rmbase( EEG, [0 50] ,[]); % 50 ms

        fprintf('\n****\nSave epoched data subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));
        SaveName = [file_type{fileno} num2str(subj_list(subj_i)) '-' num2str(sess_i) '_Preproc_ICAfull_epoched.set'];
        EEG      = pop_saveset( EEG, 'filename',SaveName,'filepath', Path2EEGsets );

        clear EEG
        ALLEEG(1:end) = [];

    end
end




%% Clean the data

% initiate or load matrix to save no. of rejected epochs and rejected components per subject and interpolated and noisy channels:
rej_epocs     = [subj_list'  nan(length(subj_list),length(sessions)*2)];
ICAcomps      = cell(length(subj_list),3);
ICAcomps(:,1) = num2cell(subj_list');
bdchns        = cell(length(subj_list),3);
bdchns(:,1)   = num2cell(subj_list');
% intrp_chans      = cell(length(subj_list),3);
% intrp_chans(:,1) = num2cell(subj_list');

rej_epocs     = table2array( readtable( [Path2EEGsets '/Overview_rejected_epochs_' char(datetime('today')) '.txt'] ) );
ICAcomps      = table2cell(  readtable( [Path2EEGsets '/Overview_ICAcomps_'        char(datetime('today')) '.txt'] ) );
bdchns        = table2cell(  readtable( [Path2EEGsets '/Overview_badchannels_'     char(datetime('today')) '.txt'] ,'Format','auto') );
% intrp_chans   = table2cell(  readtable( [Path2EEGsets '/Overview_interpolated_channels' char(datetime('today')) '.txt'] ) ); %char(datetime('yesterday')) '.txt'] ) ); %

fileno = 4;

% Loop over files
for subj_i = 6:length(subj_list)
    for sess_i = 1:length(sessions)

        fprintf('\n****\nLoad subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));
        fileName = [file_type{fileno} num2str(subj_list(subj_i)) '-' num2str(sess_i) '_Preproc_ICAfull_epoched.set'];

        % Load EEG set
        EEG      = pop_loadset('filename', fileName , 'filepath', Path2EEGsets);

        % Remove ICA components
        ALLEEG = EEG;
        CURRENTSET = 1;
        pop_eegplot( EEG, 0, 1, 1); % IC's in time domain
        EEG = pop_chanedit(EEG, 'lookup','/Users/fsmits2/Downloads/eeglab2021.0/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp');
        pop_selectcomps(EEG, 1:15 ); % IC's on topomaps and ERPs (plot first 12, only the first 20  explain a relevant amount of variance):
        % Choose which components to remove
        m = "no";
        while m ~= "yes"
            m = input('Are you ready to input components to remove? Type [yes] ','s');
        end
        m0 = -1;
        while m0 == -1
            m0 = str2double( input('How many components to remove? ','s') );
        end
        comp = [];
        if m0 > 0
            for compi = 1:m0
                comp(compi) = str2num( input(['Which component to remove? nr: ' num2str(compi) ' ' ],'s') ) ;
            end
        end
        if isempty(comp)
            ICAcomps{subj_i,sess_i+1} = '0';
        else
            ICAcomps{subj_i,sess_i+1} = comp;
            EEG = pop_subcomp( EEG, comp , 0);
        end

        % Find channels with high standard deviation to detect noisy channels
        pop_eegplot( EEG, 1, 1, 1);
        chansd     = std(EEG.data(1:30, :)');
        sortsd     = sort(chansd);
        badsd      = find( chansd > 4*mean( sortsd(1:15) ) );
        badsdchans = string({EEG.chanlocs(badsd).labels});
        if isempty(badsdchans)
            fprintf('No channels >4SD deviation\n');
        else
            fprintf('Channel with >4SD deviation: %s\n**\n', badsdchans);
        end
        % Check which channels remain bad and should not be considered further
        m3a = "no";
        while m3a ~= "yes"
            m3a = input('Ready to input channels to leave out? Type [yes] ','s');
        end
        m3 = -1;
        while m3 == -1
            m3 = str2double( input('How many channels to leave out? ','s') );
        end
        badchannels   = {[]};
        noisychannels = [];
        if m3 > 0
            for bchni = 1:m3
                badchannels{bchni}   = input(['Which channel to leave out? nr: ' num2str(bchni) ' ' ],'s') ;
                noisychannels(bchni) = find( strcmpi( badchannels{bchni}, {EEG.chanlocs.labels} ));
            end
            bdchns{subj_i,sess_i+1}  = string(badchannels); 
        else
            bdchns{subj_i,sess_i+1}  = '0';
        end
        EEG.eventdescription         = { {'Too much noise in channels: '} badchannels };

        % Leave noisy channels (previously detected and saved) out of consideration for epoch rejection
        badchannels = bdchns{subj_i, sess_i+1};
        chanarray   = 1:length(EEG.chanlocs)-2; %minus last two channels (EMG chans: HEOG, VEOG)
        if sum( strcmpi( badchannels, '0') ) < 1
            badchans = regexp(badchannels, ',', 'split');
            noisychannels = [];
            for bchni = 1:length(badchans)
                noisychannels(bchni) = find( strcmpi( badchans{bchni}, {EEG.chanlocs.labels} ));
            end
            chanarray(noisychannels) = []; % remove the noisy channels from chanarray
        end

        % Semi-automatic artifact rejection
        %     Gradient:  Specifies that the absolute difference between two adjacent sample points of data must not exceed a value (artifact of weird spikes). Starting values from Boost tutorial: Gradient: 75 μV
        %     Amplitude: Specifies that the voltage must not  exceed a certain value (artifacts like eye blinks). Starting values from Boost tutorial: Max-Min: 150 μV/200 ms
        %     Diff max-Min: Sets the threshold for the difference between the minimum and maximum voltages within the entire segment (voltage drifts). Starting values from Boost tutorial: Amplitude: -100 μV, +100 μV"
 %       ALLEEG  = EEG; CURRENTSET = 1; % Define ALLEEG and CURRENSET to enable trial rejection via pop_eegplot()
        winpnts = round(200/(1000/EEG.srate)); % points for window of 200 ms segments in the epoch
        winidx  = 1:winpnts:EEG.pnts; 
        windiff = nan(1,length(winidx));
        EEG.reject.rejmanual = zeros(1, EEG.trials); % Initialize the array for marked trials
        EEG.reject.rejmanualE = zeros(length(EEG.chanlocs), EEG.trials);

        % Loop over channels and epochs
        for ichan = chanarray(1:end-2) % exclude last two channels (VEOG & HEOG)
            for itrial = 1:EEG.trials
                gradient = max( abs( diff(EEG.data(ichan, :, itrial)) ) );
                ampliMax = max(EEG.data(ichan, :, itrial));
                ampliMin = min(EEG.data(ichan, :, itrial));
                for iwin = 1:length(winidx)-1
                    [winmin, winmax] = bounds(EEG.data( ichan, winidx(iwin):winidx(iwin)+winpnts-1, itrial));
                    windiff(iwin) = diff([winmin, winmax]);
                end
                diffV = max(windiff);

                if gradient > 50 || ampliMax > 75 || ampliMin < -75 || diffV > 100  
                    EEG.reject.rejmanual(1,itrial) = 1; % Mark the trial when a criterium is met
                    EEG.reject.rejmanualE(ichan,itrial) = 1;
                end
            end
        end

        % View the marked trials in plot
        %   Scale value to 100 and 29 epochs per window. Pay attention to VEOG.
        %   Reject the epoch around tACS artifact 
        rej_epocs(subj_i,1+sess_i) = EEG.trials; % Save total number of epochs
        find(EEG.reject.rejmanual > 0) % See the marked epoch numbers
        EEG = eeg_checkset( EEG );
        pop_eegplot( EEG, 1, 1, 0); % Plot data with marked epochs but do not immediately reject, only mark as noisy

        m0 = -1;
        while m0 == -1
            m0 = input('Ready to reject epocs? ','s');
            while isempty(m0)
                m0 = input('Ready to reject epocs? [yes]: ','s');
            end
        end

        noisyepocs = find(EEG.reject.rejmanual > 0) % see & save final series of marked epochs
        length(noisyepocs)

        m1 = [] ;
        while isempty(m1)
            m1 = input('How many epochs rejected? [enter number]: ','s');
        end

        rej_epocs(subj_i,3+sess_i) = str2double(m1);
        EEG.epochdescription       = [m1 '/' num2str(rej_epocs(subj_i,1+sess_i)) ' trials rejected'];
        EEG                        = pop_rejepoch( EEG, noisyepocs , 1);
        [ALLEEG EEG CURRENTSET]    = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off');
        EEG = eeg_checkset( EEG );
        [ALLEEG, EEG, CURRENTSET]  = eeg_store( ALLEEG, EEG );

        % Save
        fprintf('\n****\nSave clean data subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));
        SaveName = [file_type{fileno} num2str(subj_list(subj_i)) '-' num2str(sess_i) '_CleanEEG.set'];
        EEG      = pop_saveset( EEG, 'filename',SaveName,'filepath', Path2EEGsets );
        cd(path2save);
        writecell(  ICAcomps,   [Path2EEGsets '/Overview_ICAcomps_'        char(datetime('today')) '.txt'], 'Delimiter',',');
        writecell(  bdchns,     [Path2EEGsets '/Overview_badchannels_'     char(datetime('today')) '.txt'], 'Delimiter',',');
        writematrix(rej_epocs,  [Path2EEGsets '/Overview_rejected_epochs_' char(datetime('today')) '.txt'], 'Delimiter',',');

        m2 = 0;
        while m2 == 0
            m2 = input('Continue? [Y/N] ','s');
            if m2 == 'Y'
                continue
            else
                return
            end
        end

        clear EEG
        close all

        % open EEGlab
        cd('/Users/fsmits2/Downloads/eeglab2022.1')
        [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

        cd(path2scripts); % return to PITA analysis folder

    end
end




%         fprintf('\n****\nLoad subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));
%         fileName = [file_type{fileno} num2str(subj_list(subj_i)) '-' num2str(sess_i) '_CleanEEG.set'];
% 
%         % Load EEG set
%         EEG      = pop_loadset('filename', fileName , 'filepath', Path2EEGsets);
% 
%         % Reverse eyes-open eyes-closed codes for subj 262 session 1
%         if subj_list(subj_i)==262 && sess_i==1
%             for trigi = 1:length(EEG.event)
%                 if EEG.event(trigi).type=='11'
%                     EEG = pop_editeventvals(EEG,'changefield', {trigi 'type' '0'});
%                 elseif EEG.event(trigi).type=='22'
%                     EEG = pop_editeventvals(EEG,'changefield', {trigi 'type' '11'});
%                 end
%             end
%             for trigi = 1:length(EEG.event)
%                 if EEG.event(trigi).type=='0'
%                     EEG = pop_editeventvals(EEG,'changefield', {trigi 'type' '22'});
%                 else
%                     continue
%                 end
%             end
%         end








% %% Re-do cleaning for correctly rejecting epocs and saving set with clean epochs only
% 
% % initiate or load matrix to save no. of rejected epochs and rejected components per subject and interpolated and noisy channels:
% rej_epocs     = [subj_list'  nan(length(subj_list),length(sessions)*2)];
% 
% ICAcomps      = table2cell(  readtable( [Path2EEGsets '/Overview_ICAcomps_20-Jan-2023.txt'] ) );
% bdchns        = table2cell(  readtable( [Path2EEGsets '/Overview_badchannels_20-Jan-2023.txt'] ,'Format','auto') );
% 
% rej_epocs     = table2array( readtable( [Path2EEGsets '/Overview_rejected_epochs_' char(datetime('today')) '.txt'] ) );
% 
% fileno = 4;
% 
% % Loop over files
% for subj_i = 31:length(subj_list)
%     for sess_i = 1:length(sessions)
% 
%         fprintf('\n****\nLoad subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));
%         fileName = [file_type{fileno} num2str(subj_list(subj_i)) '-' num2str(sess_i) '_CleanEEG.set'];
% 
%         % Load EEG set
%         EEG      = pop_loadset('filename', fileName, 'filepath', Path2EEGsets);
% 
%         [ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );
%         EEG = eeg_checkset( EEG );
% 
% %         % Leave noisy channels (previously detected and saved) out of consideration for epoch rejection
% %         badchannels = bdchns{subj_i, sess_i+1};
% %         chanarray   = 1:length(EEG.chanlocs)-2; %minus last two channels (EMG chans: HEOG, VEOG)
% %         if sum( strcmpi( badchannels, '0') ) < 1
% %             badchans = regexp(badchannels, ',', 'split');
% %             noisychannels = [];
% %             for bchni = 1:length(badchans)
% %                 noisychannels(bchni) = find( strcmpi( badchans{bchni}, {EEG.chanlocs.labels} ));
% %             end
% %             chanarray(noisychannels) = []; % remove the noisy channels from chanarray
% %         end
% % 
% %         % Semi-automatic artifact rejection
% %         %     Gradient:  Specifies that the absolute difference between two adjacent sample points of data must not exceed a value (artifact of weird spikes). Starting values from Boost tutorial: Gradient: 75 μV
% %         %     Amplitude: Specifies that the voltage must not  exceed a certain value (artifacts like eye blinks). Starting values from Boost tutorial: Max-Min: 150 μV/200 ms
% %         %     Diff max-Min: Sets the threshold for the difference between the minimum and maximum voltages within the entire segment (voltage drifts). Starting values from Boost tutorial: Amplitude: -100 μV, +100 μV"
% %  %       ALLEEG  = EEG; CURRENTSET = 1; % Define ALLEEG and CURRENSET to enable trial rejection via pop_eegplot()
% %         winpnts = round(200/(1000/EEG.srate)); % points for window of 200 ms segments in the epoch
% %         winidx  = 1:winpnts:EEG.pnts; 
% %         windiff = nan(1,length(winidx));
% %         EEG.reject.rejmanual = zeros(1, EEG.trials); % Initialize the array for marked trials
% %         EEG.reject.rejmanualE = zeros(length(EEG.chanlocs), EEG.trials);
% % 
% %         % Loop over channels and epochs
% %         for ichan = chanarray(1:end-2) % exclude last two channels (VEOG & HEOG)
% %             for itrial = 1:EEG.trials
% %                 gradient = max( abs( diff(EEG.data(ichan, :, itrial)) ) );
% %                 ampliMax = max(EEG.data(ichan, :, itrial));
% %                 ampliMin = min(EEG.data(ichan, :, itrial));
% %                 for iwin = 1:length(winidx)-1
% %                     [winmin, winmax] = bounds(EEG.data( ichan, winidx(iwin):winidx(iwin)+winpnts-1, itrial));
% %                     windiff(iwin) = diff([winmin, winmax]);
% %                 end
% %                 diffV = max(windiff);
% % 
% %                 if gradient > 50 || ampliMax > 75 || ampliMin < -75 || diffV > 100  
% %                     EEG.reject.rejmanual(1,itrial) = 1; % Mark the trial when a criterium is met
% %                     EEG.reject.rejmanualE(ichan,itrial) = 1;
% %                 end
% %             end
% %         end
% % 
% %         [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
% % 
% %         % View the marked trials in plot
% %         %   Scale value to 100 and 29 epochs per window. Pay attention to VEOG.
% %         %   Reject the epoch around tACS artifact 
% %         rej_epocs(subj_i,1+sess_i) = EEG.trials; % Save total number of epochs
% %         find(EEG.reject.rejmanual > 0) % See the marked epoch numbers
% %         EEG = eeg_checkset( EEG );
% %         pop_eegplot( EEG, 1, 1, 0); % Plot data with marked epochs but do not immediately reject, only mark as noisy
% % 
% % 
% %         m0 = -1;
% %         while m0 == -1
% %             m0 = input('Ready to reject epocs? ','s');
% %             while isempty(m0)
% %                 m0 = input('Ready to reject epocs? [yes]: ','s');
% %             end
% %         end
% % 
% %         noisyepocs = find(EEG.reject.rejmanual > 0) % see & save final series of marked epochs
% % 
% %         m1 = [] ;
% %         while isempty(m1)
% %             m1 = input('How many epochs rejected? [enter number]: ','s');
% %         end
% % 
% %         rej_epocs(subj_i,3+sess_i) = str2double(m1);
% %         EEG.epochdescription       = [m1 '/' num2str(rej_epocs(subj_i,1+sess_i)) ' trials rejected'];
% %         EEG                        = pop_rejepoch( EEG, noisyepocs , 1);
% %         [ALLEEG EEG CURRENTSET]    = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off');
% %         EEG = eeg_checkset( EEG );
% %         [ALLEEG, EEG, CURRENTSET]  = eeg_store( ALLEEG, EEG );
% 
% 
%         % Save
%         fprintf('\n****\nSave clean data subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));
%         SaveName = [file_type{fileno} num2str(subj_list(subj_i)) '-' num2str(sess_i) '_CleanEEG_CleanEpochs.set'];
%         EEG      = pop_saveset( EEG, 'filename',SaveName,'filepath', Path2EEGsets );
%         writematrix(rej_epocs,  [Path2EEGsets '/Overview_rejected_epochs_' char(datetime('today')) '.txt'], 'Delimiter',',');
%        
%         m2 = 0;
%         while m2 == 0
%             m2 = input('Continue? [Y/N] ','s');
%             if m2 == 'Y'
%                 continue
%             else
%                 return
%             end
%         end
% 
%         clear EEG
%         close all
% 
%         % open EEGlab
%         cd('/Users/fsmits2/Downloads/eeglab2022.1')
%         [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
% 
%         cd(path2scripts); % return to PITA analysis folder
% 
%     end
% end
% 

