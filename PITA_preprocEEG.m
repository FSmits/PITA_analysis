%% EEG pre-processing - PITA - all data

clear
close all

%% open EEGlab

cd('/Users/fsmits2/Downloads/eeglab2021.0')
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% return to PITA analysis folder
cd('/Users/fsmits2/Documents/PITA_analysis');


%% set paths to data

Path2EEGbdf  = '/Users/fsmits2/Downloads/1 EEG data tacs-eeg/';
Path2EEGsets = '/Users/fsmits2/Downloads/1 EEG data tacs-eeg/3 EEG sets processed';


%% enter subject identification numbers

subj_list =[669	557 363	638	602	989	383	502	733	442	575	913 710	262 ...
    752 227	565	626 334	362	600	121	319 923	915	298	202	692	275	...
    508 291	803	755	681	876	134	559	396	818	601	297	524	883	193	642];

% NOTE: 710-1 NO STIMULATION! (but =sham session)

sessions  =  [1 2];


%% Original trigger codes

StartRec        = 254;
StopRec         = 255;
start_tACS      = 5; %(2:20) identical to stop_EEG minus(1:19)
stop_tACS       = 6; %(1:20) identical to start_EEG (1:20)
start_EEG       = 7;
stop_EEG        = 8;
post_tACS_EEG   = 9;


%% enter filenames

recording1 = 'restingstate-pretACS-';
recording2 = 'Encoding-';
recording3 = 'TACSEEG-';
recording4 = 'restingstate-posttACS-';
recording5 = 'Retrieval-';

file_type = {recording1, recording2, recording3, recording4, recording5};

%% read EEG file

for subj_i = 1:length(subj_list)

    for sess_i = 1:length(sessions)

        fprintf('\n****\nStart processing subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));

        fileName = fullfile(Path2EEGbdf, [file_type{3} num2str(subj_list(subj_i)) '-' num2str(sess_i) '.bdf']);

        % Load raw bdf data file via EEGlab
        EEG = pop_biosig( fileName );


        % Enter data to the EEG structure
        EEG.filename = fileName;
        EEG.setname  = fileName;
        EEG.subject  = subj_list(subj_i);
        EEG.session  = sess_i;


        % Remove non-recorded channels
        % F3 F4 (tACS electrode locations) and EXG7 EXG8
        EEG = pop_select(EEG, 'nochannel', {'F3', 'F4', 'EXG7', 'EXG8'});


        % Re-code events
        % [Why? BioSemi/Computer settings resulted in changes in the recorded trigger codes relative to the originally programmed triggers codes. These changes are unfortunately not exactly the same across subjects.]

        % remove added trigger text like 'condition' and 'artifact'
        for ev_i = 1:length({EEG.event.type})
            EEG.event(ev_i).type = strrep( EEG.event(ev_i).type, 'condition ', '' );
            EEG.event(ev_i).type = strrep( EEG.event(ev_i).type, 'artifact', '' );
        end

        % remove trigger 256
        trig_256 = find(strcmpi( {EEG.event.type}, '256' ));
        EEG = pop_editeventvals(EEG,'delete', trig_256);

        % find or create trigger for tACS onset
        % (original code = '5'. Other possibile codes = '8' , ':Failing electrode' , ':breathing'.)
        if any( subj_list(subj_i) == [638 275 262] )  &&  sess_i==1
            EEG = pop_editeventvals(EEG,'insert',{1,[],[],[],[]},'changefield',{1,'type','5'},'changefield',{1,'edftype','5'},'changefield',{1,'latency',0.01});
        end
        tACS_begin = sort([find(strcmpi( {EEG.event.type}, '5' )), ...
            find(strcmpi( {EEG.event.type}, '8' )),...
            find(strcmpi( {EEG.event.type}, ':Failing electrode' )),...
            find(strcmpi( {EEG.event.type}, ':breathing' )) ]);

        EEG_end     = tACS_begin(end);

        % exceptions
        if subj_list(subj_i) == 669  &&  sess_i == 2
            tACS_begin([2 12 14 24])  = []; % events 4, 31, 36, 64
        elseif subj_list(subj_i) == 915  &&  sess_i == 1
            tACS_begin([2 12 14 24])  = []; % events 5 32 37 65
        elseif subj_list(subj_i) == 202  &&  sess_i == 1
            tACS_begin([2 12 14 24])  = []; % events 4 33 38 66
        elseif subj_list(subj_i) == 298  &&  sess_i == 1
            tACS_begin([2 4 14 24])  = []; % events 4 9 38 67
        elseif subj_list(subj_i) == 638  &&  sess_i == 2
            tACS_begin([15 17 23])    = []; % events 44, 49, 65
        elseif subj_list(subj_i) == 710  &&  sess_i == 1
            tACS_begin([7 9 23])  = []; % events 19, 24,65
        elseif subj_list(subj_i) == 227  &&  sess_i == 2
            tACS_begin([5 15 17 24])  = []; % events 14 42 47 66
        elseif subj_list(subj_i) == 334  &&  sess_i == 1
            tACS_begin([3 5 17 24])  = []; % events 7 12 48 68
        elseif subj_list(subj_i) == 319  &&  sess_i == 2
            tACS_begin([1 11 13 15 25])  = []; % events 1 30 35 39 67
        elseif subj_list(subj_i) == 923  &&  sess_i == 1
            tACS_begin([6 15 17 19 25])  = []; % events 16 42 46 51 67
        elseif subj_list(subj_i) == 298  &&  sess_i == 2
            tACS_begin([9 11 21 23 25])  = []; % events 25 30 57 62 66
        elseif subj_list(subj_i) == 275  &&  sess_i == 2
            tACS_begin([4 6 23])  = []; % event 10 15 64
        elseif subj_list(subj_i) == 508  &&  sess_i == 1
            tACS_begin([16 18 23 ])  = []; % event 47 51 64
        elseif subj_list(subj_i) == 291  &&  sess_i == 1
            tACS_begin([3 13 24 25])  = []; % event 6 35 68 69
        elseif subj_list(subj_i) == 876  &&  sess_i == 2
            tACS_begin([2 22])  = []; % event 5 63
        elseif subj_list(subj_i) == 559  &&  sess_i == 1
            tACS_begin([16 18 23])  = []; % event 46 51 64
        elseif subj_list(subj_i) == 396  &&  sess_i == 1
            tACS_begin([6 16 23])  = []; % event 16 45 65
        end

        if length(tACS_begin)>19
            tACS_begin  = tACS_begin(1:20);
        end

        % find trigger for tACS offset
        % (original code = '6'. Other possibile codes = '7' , ':sweat' , ':50/60 Hz mains interference'.)
        tACS_end = sort([find(strcmpi( {EEG.event.type}, '6' )), ...
            find(strcmpi( {EEG.event.type}, '7' )), ...
            find(strcmpi( {EEG.event.type}, ':50/60 Hz mains interference' )), ...
            find(strcmpi( {EEG.event.type}, ':sweat' )) ]);
        if length(tACS_end) > 30
            tACS_end = tACS_end(2:2:40);
        end

        % re-code to original code for tACS onsets and offsets
        for trig_i = 1:length(tACS_end)
            EEG = pop_editeventvals(EEG,'changefield', {tACS_begin(trig_i) 'type' '5'});
            EEG = pop_editeventvals(EEG,'changefield', {tACS_end(trig_i) 'type' '6'});
        end


        % Cut out data during stimulation (tACS)
        for ni = flip(1:length(tACS_begin))
            EEG = eeg_eegrej(EEG, [EEG.event(tACS_begin(ni)).latency,  EEG.event(tACS_end(ni)).latency]); % remove data during tACS
            [ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG); % store changes
        end


        % Cut out data before first tACS-EEG & after last tACS-EEG period
        EEG_begin = find(strcmpi( {EEG.event.type}, 'boundary' ));
        EEG_end   = sort([find(strcmpi( {EEG.event.type}, '8' )), find(strcmpi( {EEG.event.type}, ':breathing' ))]);
        EEG_end   = EEG_end(1,end); % in case multiple trigger codes '8' are available, pick last one

        EEG       = eeg_eegrej(EEG, [EEG.event(EEG_end).latency + 0.05,   EEG.pnts(end)]);
        EEG       = eeg_eegrej(EEG, [0,   EEG.event(EEG_begin(1)).latency - 0.05]);


        %         % Inspect data
        %         pop_eegplot( EEG, 1, 1, 1);


        % Save
        fprintf('\n****\nSave processed subject %i session %i\n****\n\n', ...
            subj_list(subj_i), sessions(sess_i));
        SaveName = [file_type{3} num2str(subj_list(subj_i)) '-' num2str(sess_i) '_RawEEG.set'];
        EEG = pop_saveset( EEG, 'filename',SaveName,'filepath', Path2EEGsets );

        clear EEG
        ALLEEG(1:end) = [];

    end

end


%% Clean each individual dataset

for subj_i = 1:2 %length(subj_list)

   % for sess_i = 1:length(sessions)

        fprintf('\n****\nStart processing subject %i session %i\n****\n\n', subj_list(subj_i), sessions(sess_i));

        fileName = ['TACSEEG-' num2str(subj_list(subj_i)) '-' num2str(sess_i)];

        % Load EEG set
        EEG      = pop_loadset('filename', [fileName, '_RawEEG.set'], 'filepath', Path2EEGsets);

        % Plot data and check if tACS artifact is visible. If not - check if stimulation was carried out correctly.
        pop_eegplot( EEG, 1, 1, 1);

        m = 0 ;
        while m == 0
            m = input('Do you see tACS artifacts? Y/N:','s');
            if m == 'N'
                writecell( { ['tACS-EEG-' num2str(subj_list(subj_i)) '-' num2str(sess_i)] } , ...
                    ['/Users/fsmits2/Downloads/1 EEG data tacs-eeg/3 EEG sets processed/' ...
                    ['noStim-' num2str(subj_list(subj_i)) '-' num2str(sess_i) '.txt']]);
                continue
            elseif m == 'Y'


                hansie = 199



            end
        end


   % end
end





