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
cd('/Users/fsmits2/Downloads/eeglab2021.0')
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

cd('/Users/fsmits2/Documents/PITA_analysis'); % return to PITA analysis folder

% set paths to data
Path2EEGbdf  = '/Users/fsmits2/Downloads/1 EEG data tacs-eeg/';
Path2EEGsets = '/Users/fsmits2/Downloads/1 EEG data tacs-eeg/3 EEG sets processed';

% enter subject names
subj_list =[669	557 363	638	989	383	502	733	442	575	710	262 ...
    752 227	565	362	600	121	319 923	915	298	202	692	275	...
    508 291	803	755	681	876	134	559	818	601	524	883	193	642];
sessions  = [1 2];


% enter filenames
recording1 = 'restingstate-pretACS-';
recording2 = 'Encoding-';
recording3 = 'TACSEEG-';
recording4 = 'restingstate-posttACS-';
recording5 = 'Retrieval-';
file_type = {recording1, recording2, recording3, recording4, recording5};


%% read EEG file

% Original trigger codes
StartRec        = 254;
StopRec         = 255;
start_tACS      = 5; %(2:20) identical to stop_EEG minus(1:19)
stop_tACS       = 6; %(1:20) identical to start_EEG (1:20)
start_EEG       = 7;
stop_EEG        = 8;
post_tACS_EEG   = 9;

% Initiate time/period-related triggers
secs      = 30; % 30 seconden data na elke tACS stimulatie
stims     = 20; % 20 tACS stimulaties in totaal
trig_base = repmat(0.01:0.01:0.30,[stims,1]);
stim_mat  = repmat(1:stims,[secs,1])';
trigs     = trig_base + stim_mat; % Trigger names are: #stimulation as integer, #second of data following that stimulation as decimal

% Loop over files
for subj_i = 1:length(subj_list)
    for sess_i = 1:length(sessions)

        fprintf('\n****\nStart processing subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));

        fileName = fullfile(Path2EEGbdf, [file_type{3} num2str(subj_list(subj_i)) '-' num2str(sess_i) '.bdf']);

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

        % find or create trigger for tACS onset   (original code = '5'. Other possibile codes = '8' , ':Failing electrode' , ':breathing'.)
        if any( subj_list(subj_i) == [638 275 262] )  &&  sess_i==1
            EEG = pop_editeventvals(EEG,'insert',{1,[],[],[],[]},'changefield',{1,'type','5'},'changefield',{1,'edftype','5'},'changefield',{1,'latency',0.01});
        end
        tACS_begin = sort([find(strcmpi( {EEG.event.type}, '5' )), ...
            find(strcmpi( {EEG.event.type}, '8' )),...
            find(strcmpi( {EEG.event.type}, ':Failing electrode' )),...
            find(strcmpi( {EEG.event.type}, ':breathing' )) ]);

        % find trigger for tACS offset   (original code = '6'. Other possibile codes = '7' , ':sweat' , ':50/60 Hz mains interference'.)
        tACS_end = sort([find(strcmpi( {EEG.event.type}, '6' )), ...
            find(strcmpi( {EEG.event.type}, '7' )), ...
            find(strcmpi( {EEG.event.type}, ':50/60 Hz mains interference' )), ...
            find(strcmpi( {EEG.event.type}, ':sweat' )) ]);

        % when too many trigger codes are found, select codes that show >90-second trigger-to-trigger delay (tACS_begin and tACS_end triggers are separated by ~92 seconds: 60 seconds stimulation + 30 second EEG recording period + 2 seconds ramp-up/ramp-down stimulation)
        if length(tACS_begin) > 21
            evlat90_begin = [];
            evlat = [];
            for eventi = 1:length(tACS_begin)
                evlat(eventi) = EEG.event(tACS_begin(eventi)).latency / EEG.srate; %save event latencies in array
            end
            evlat90_begin  = find( diff(evlat) > 90) + 1; %find which events have a latency difference of >90 seconds
            tACS_begin     = [tACS_begin(1) tACS_begin(evlat90_begin)];
        end
        if length(tACS_end) > 21
            evlat90_end = [];
            evlat = [];
            for eventi = 1:length(tACS_end)
                evlat(eventi) = EEG.event(tACS_end(eventi)).latency / EEG.srate; %save event latencies in array
            end
            evlat90_end  = find( diff(evlat) > 90) + 1; %find which events have a latency difference of >90 seconds
            tACS_end     = [tACS_end(1) tACS_end(evlat90_end)];
        end

        % re-code to original code for tACS onsets and offsets
        for trig_i = 1:length(tACS_end)
            EEG = pop_editeventvals(EEG,'changefield', {tACS_begin(trig_i) 'type' 'tACS_start'});
            EEG = pop_editeventvals(EEG,'changefield', {tACS_end(trig_i)   'type' 'tACS_stop'});
        end

        % insert time-related indices as event codes
        bndrs     = find(strcmpi( {EEG.event.type}, 'tACS_stop' ));
        bndrs_lat = [EEG.event(bndrs).latency];
        % insert indexing trigger code every full second + little delay of 5 ms so that the indexing triggers do not overlap with 1-second epoch cuts
        for i_EEGs = 1:length(bndrs)
            for i_secs = 1:secs
                EEG = pop_editeventvals(EEG, 'insert',{bndrs(i_EEGs)+i_secs,[],[],[],[]},...
                    'changefield',{bndrs(i_EEGs)+i_secs ,'type',    trigs(i_EEGs,i_secs) },...
                    'changefield',{bndrs(i_EEGs)+i_secs ,'edftype', trigs(i_EEGs,i_secs) },...
                    'changefield',{bndrs(i_EEGs)+i_secs ,'latency', bndrs_lat(i_EEGs)/EEG.srate + i_secs}); % for latency event, latencies are in millisecond compared to the time locking event, not in data samples.
            end
        end

        % -- Cut out data during stimulation (tACS)
        tACS_begin = sort(find(strcmpi( {EEG.event.type}, 'tACS_start' )));
        tACS_end   = sort(find(strcmpi( {EEG.event.type}, 'tACS_stop' )));
        for ni = flip(1:length(tACS_begin))
            EEG = eeg_eegrej(EEG, [EEG.event(tACS_begin(ni)).latency,  EEG.event(tACS_end(ni)).latency]); % remove data during tACS
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG); % store changes
        end

        % -- Cut out data before first tACS-EEG & after last tACS-EEG period
        EEG_begin = find(strcmpi( {EEG.event.type}, 'boundary' ));
        EEG_end   = sort([find(strcmpi( {EEG.event.type}, '8' )), find(strcmpi( {EEG.event.type}, ':breathing' ))]);
        EEG_end   = EEG_end(1,end); % in case multiple trigger codes '8' are available, pick last one
        EEG       = eeg_eegrej(EEG, [EEG.event(EEG_end).latency + 0.05,   EEG.pnts(end)]);
        EEG       = eeg_eegrej(EEG, [0,   EEG.event(EEG_begin(1)).latency - 0.05]);

        %         pop_eegplot( EEG, 1, 1, 1); % Inspect data

        % -- Save
        fprintf('\n****\nSave processed subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));
        SaveName = [file_type{3} num2str(subj_list(subj_i)) '-' num2str(sess_i) '_RawEEG.set'];
        EEG = pop_saveset( EEG, 'filename',SaveName,'filepath', Path2EEGsets );

        clear EEG
        ALLEEG(1:end) = [];
    end
end



%% Pre-processing steps: re-reference, downsample, filter, create bipolar EOG channels

for subj_i = 1:length(subj_list)
    for sess_i = 1:length(sessions)

        fprintf('\n****\nStart pre-processing subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));

        fileName = ['TACSEEG-' num2str(subj_list(subj_i)) '-' num2str(sess_i)];

        % Load EEG set
        EEG      = pop_loadset('filename', [fileName, '_RawEEG.set'], 'filepath', Path2EEGsets);

        % Check if all EEG-recording episodes are present
        m0 = -1 ;
        fprintf([ num2str(subj_list(subj_i)) ' session ' num2str(sessions(sess_i)) ' -- Number of events: ' num2str(length(EEG.event)) ' \n '])
        while m0 == -1
            if length(EEG.event) < 600
                m0 = input('Continue? Y/N: ','s');
                if m0 == 'Y'
                    continue
                else 
                    m0 = 0;
                    return
                end
            else
                m0 = 0;
                continue
            end
        end

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

        %         pop_eegplot( EEG, 1, 1, 1); % Inspect data

        % Save
        fprintf('\n****\nSave pre-processed subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));
        SaveName = [file_type{3} num2str(subj_list(subj_i)) '-' num2str(sess_i) '_PreprocEEG.set'];
        EEG = pop_saveset( EEG, 'filename',SaveName,'filepath', Path2EEGsets );

        clear EEG
        ALLEEG(1:end) = [];
    end
end



%% Interpolate channels and detect gel bridges

% initiate or load matrix/cellarray for gel bridges and interpolated channels per subject:
gelbridges       = [subj_list'  nan(length(subj_list),length(sessions))];
intrp_chans      = cell(length(subj_list),3);
intrp_chans(:,1) = num2cell(subj_list');

gelbridges       = table2array( readtable( [Path2EEGsets '/Overview_gelbridges_'           char(datetime('today')) '.txt'] ) ); %char(datetime('yesterday')) '.txt'] ) ); %
intrp_chans      = table2cell(  readtable( [Path2EEGsets '/Overview_interpolated_channels' char(datetime('today')) '.txt'] ) ); %char(datetime('yesterday')) '.txt'] ) ); %

% Loop over files
for subj_i = 1:length(subj_list)
    for sess_i = 1:length(sessions)

        fprintf('\n****\nLoad subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));
        fileName = ['TACSEEG-' num2str(subj_list(subj_i)) '-' num2str(sess_i)];

        % Load and plot EEG set
        EEG      = pop_loadset('filename', [fileName, '_PreprocEEG.set'], 'filepath', Path2EEGsets);
        pop_eegplot(EEG,1,1,1);

        % Interpolate bad channels
        k0 = [];
        m0 = -1;
        while m0 == -1
            m0 = str2double( input('How many bad channels (interpolation needed?) ','s') );
        end
        m0s = num2cell(nan(1,m0));
        k0 = input('Gel bridge? [1 = yes / 0 = no] ');
        if k0 == 1
            gelbridges(subj_i,sess_i+1) = 1;
        else
            gelbridges(subj_i,sess_i+1) = 0;
        end
        badchan = {[]};
        if m0 > 0
            for badchani = 1:m0
                badchan{badchani} = input(['Which channel to interpolate? nr ' num2str(badchani) ' ' ],'s');
                chan2interp       = find( strcmpi( badchan{badchani}, {EEG.chanlocs.labels} ));
                % Enter channel locations
                EEG = pop_chanedit(EEG, 'lookup','/Users/fsmits2/Downloads/eeglab2021.0/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp');
                EEG = pop_interp(EEG, chan2interp, 'spherical');
            end
        end
        if isempty(cell2mat(badchan))
            intrp_chans{subj_i,sess_i+1} = '0';
        else
            intrp_chans{subj_i,sess_i+1} = string(badchan);
        end
        EEG.eventdescription         = { {'Interpolated channels: '} badchan };
       
        % Save
        fprintf('\n****\nSave interpolated channels data subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));
        SaveName = [file_type{3} num2str(subj_list(subj_i)) '-' num2str(sess_i) '_intrp.set'];
        EEG      = pop_saveset( EEG, 'filename',SaveName,'filepath', Path2EEGsets );
        writematrix(gelbridges, [Path2EEGsets '/Overview_gelbridges_'           char(datetime('today')) '.txt'], 'Delimiter',';'); %char(datetime('yesterday')) '.txt'], 'Delimiter',';'); %
        writecell(intrp_chans,  [Path2EEGsets '/Overview_interpolated_channels' char(datetime('today')) '.txt'], 'Delimiter',';'); %char(datetime('yesterday')) '.txt'], 'Delimiter',';'); %

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
        ALLEEG(1:end) = [];
    end
end


%% ICA

% Loop over files
for subj_i = 1:length(subj_list)
    for sess_i = 1:length(sessions)

        fprintf('\n****\nLoad subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));
        fileName = ['TACSEEG-' num2str(subj_list(subj_i)) '-' num2str(sess_i)];

        % Load EEG set
        EEG      = pop_loadset('filename', [fileName, '_intrp.set'], 'filepath', Path2EEGsets);

        % Find Eyeblink components with ICA (runica):
        EEG = pop_runica(EEG, 'icatype', 'runica', 'extended',1,'interrupt','on');

        % Save
        fprintf('\n****\nSave ICA data subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));
        SaveName = [file_type{3} num2str(subj_list(subj_i)) '-' num2str(sess_i) '_intrp_ICAfull.set'];
        EEG      = pop_saveset( EEG, 'filename',SaveName,'filepath', Path2EEGsets );
        writecell(ICAcomps,  [Path2EEGsets '/Overview_ICAcomps_' char(datetime('today')) '.txt'], 'Delimiter',';');
  
        clear EEG
        ALLEEG(1:end) = [];

    end
end


%% Clean the data

% initiate or load matrix to save no. of rejected epochs and rejected components per subject:
rej_epocs     = [subj_list'  nan(length(subj_list),length(sessions)*2)];
ICAcomps      = cell(length(subj_list),3);
ICAcomps(:,1) = num2cell(subj_list');

rej_epocs     = table2array( readtable( [Path2EEGsets '/Overview_rejected_epochs_'      char(datetime('today')) '.txt'] ) );
ICAcomps      = table2cell(  readtable( [Path2EEGsets '/Overview_ICAcomps_' char(datetime('today')) '.txt'] ) );

% Loop over files
for subj_i = 23:length(subj_list)
    for sess_i = 1:length(sessions)

        fprintf('\n****\nLoad subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));
        fileName = ['TACSEEG-' num2str(subj_list(subj_i)) '-' num2str(sess_i)];

        % Load EEG set
        EEG      = pop_loadset('filename', [fileName, '_intrp_ICAfull.set'], 'filepath', Path2EEGsets);

        % Remove ICA components
        ALLEEG = EEG;
        CURRENTSET = 1;
        pop_eegplot( EEG, 0, 1, 1); % IC's in time domain
        EEG = pop_chanedit(EEG, 'lookup','/Users/fsmits2/Downloads/eeglab2021.0/plugins/dipfit/standard_BESA/standard-10-5-cap385.elp');
        pop_selectcomps(EEG, 1:12 ); % IC's on topomaps and ERPs (plot first 12, only the first 20  explain a relevant amount of variance):
        % Choose which components to remove
        m = "no";
        while m ~= "yes"
            m = input('Are you ready to input components to remove? Type [yes]','s');
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

        % Find and remove irrelevant triggers in EEG-recording periods
        trigarray = [];
        for eventi = 1:length(EEG.event)
            trigarray(eventi) = str2double(EEG.event(eventi).type);
        end
        rejects_char = find( isnan( trigarray) );
        rejects_num  = find( trigarray == 8 | trigarray == 7 | trigarray == 6 | trigarray == 5 );
        trigarray( [rejects_char rejects_num] ) = [];

        % Divide into 1-sec epochs
        EEG = pop_epoch( EEG, num2cell(trigarray),  [0  1], 'epochinfo', 'yes');
        % Remove baseline mean
        EEG = pop_rmbase( EEG, [0 50] ,[]); % 50 ms

        % Semi-automatic artifact rejection
        %     Gradient:  Specifies that the absolute difference between two adjacent sample points of data must not exceed a value (artifact of weird spikes). Starting values from Boost tutorial: Gradient: 75 μV
        %     Amplitude: Specifies that the voltage must not  exceed a certain value (artifacts like eye blinks). Starting values from Boost tutorial: Max-Min: 150 μV/200 ms
        %     Diff max-Min: Sets the threshold for the difference between the minimum and maximum voltages within the entire segment (voltage drifts). Starting values from Boost tutorial: Amplitude: -100 μV, +100 μV"
        ALLEEG  = EEG; CURRENTSET = 1; % Define ALLEEG and CURRENSET to enable trial rejection via pop_eegplot()
        winpnts = round(200/(1000/EEG.srate)); % points for window of 200 ms segments in the epoch
        winidx  = 1:winpnts:EEG.pnts; 
        windiff = nan(1,length(winidx));
        EEG.reject.rejmanual = zeros(1, EEG.trials); % Initialize the array for marked trials
        EEG.reject.rejmanualE = zeros(length(EEG.chanlocs), EEG.trials);

        % Loop over channels and epochs
        for ichan = 1:length(EEG.chanlocs)-2 % exclude last two channels (VEOG & HEOG)
            for itrial = 1:EEG.trials
                gradient = max( abs( diff(EEG.data(ichan, :, itrial)) ) );
                ampliMax = max(EEG.data(ichan, :, itrial));
                ampliMin = min(EEG.data(ichan, :, itrial));
                for iwin = 1:length(winidx)-1
                    [winmin, winmax] = bounds(EEG.data( ichan, winidx(iwin):winidx(iwin)+winpnts-1, itrial));
                    windiff(iwin) = diff([winmin, winmax]);
                end
                diffV = max(windiff);

                if gradient > 75 || ampliMax > 100 || ampliMin < -100 || diffV > 150  
                    EEG.reject.rejmanual(1,itrial) = 1; % Mark the trial when a criterium is met
                    EEG.reject.rejmanualE(ichan,itrial) = 1;
                end
            end
        end

        % View the marked trials in plot
        %   Scale value to 80 and 25 epochs per window. Pay attention to VEOG.
        %   Reject the epoch around tACS artifact 
        rej_epocs(subj_i,1+sess_i) = EEG.trials; % Save total number of epochs
        find(EEG.reject.rejmanual > 0) % See the marked epoch numbers
        pop_eegplot( EEG, 1, 1, 1); % Plot data with marked epochs
        
        m1 = -1 ;
        while m1 == -1
            m1 = input('How many epochs rejected?: ','s');
            while isempty(m1)
                m1 = input('How many epochs rejected? [enter number]: ','s');
            end
        end
        rej_epocs(subj_i,3+sess_i) = str2double(m1);
        EEG.epochdescription         = [m1 '/' num2str(rej_epocs(subj_i,1+sess_i)) ' trials rejected'];

        % Save
        fprintf('\n****\nSave clean data subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));
        SaveName = [file_type{3} num2str(subj_list(subj_i)) '-' num2str(sess_i) '_CleanEEG_inclBlinks.set'];
        EEG      = pop_saveset( EEG, 'filename',SaveName,'filepath', Path2EEGsets );
        writematrix(rej_epocs,  [Path2EEGsets '/Overview_rejected_epochs_'      char(datetime('today')) '.txt'], 'Delimiter',';');
        
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
        ALLEEG(1:end) = [];
    end
end
