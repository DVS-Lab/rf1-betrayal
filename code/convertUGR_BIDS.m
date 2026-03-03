%% convert UGR raw behavioral to BIDS events (rf1-betrayal version)

%% ---- DEFINE PROJECT ROOT ----
project_root = '/ZPOOL/data/projects/rf1-betrayal';

basedir = project_root;   % where BIDS folder will live
datadir = '/ZPOOL/data/projects/rf1-sra/stimuli/Scan-Lets_Make_A_Deal';

%% ---- LOAD SUBJECT LIST ----
sublist_file = fullfile(project_root, 'code', 'sublist_n139.txt');
assert(isfile(sublist_file), 'Subject list file not found.')

subjects = readmatrix(sublist_file);
subjects = subjects(~isnan(subjects));  % remove blank rows

fprintf('Loaded %d subjects\n', length(subjects));

%% ---- LOOP THROUGH SUBJECTS ----
for s = 1:length(subjects)

    subID = subjects(s);
    subID_str = sprintf('%05d', subID);

    for r = 0:1

        %% ---- LOAD RAW DATA ----
        indata = fullfile(datadir, 'logs', subID_str, ...
            sprintf('sub-%05d_task-ultimatum_run-%01d_raw.csv', subID, r));

        if ~exist(indata, 'file')
            fprintf('Missing: %s\n', indata);
            continue
        end

        T = readtable(indata);

        %% ---- VALIDATE ROW COUNT ----
        valid_rows = all(~ismissing(T), 2);
        num_valid_rows = sum(valid_rows);

        if num_valid_rows ~= 48
            fprintf('Invalid row count (n=%d): %s\n', num_valid_rows, indata);
            continue
        end

        %% ---- EXTRACT COLUMNS ----
        decision_onset = T.decision_onset;
        onset = T.cue_Onset;
        RT = T.rt;
        Block = T.Block;
        Endowment = T.Endowment;
        response = T.resp;
        L_Option = T.L_Option;
        R_Option = T.R_Option;

        offer = max([L_Option R_Option],[],2);

        if isa(decision_onset,'cell')
            fprintf('Cell array detected (incomplete data): %s\n', indata);
            continue
        end

        %% ---- CREATE BIDS OUTPUT FILE ----
        fname = sprintf('sub-%05d_task-ugr_run-%01d_events.tsv', subID, r+1);
        output = fullfile(basedir, 'bids', ['sub-' subID_str], 'func');

        if ~exist(output,'dir')
            mkdir(output)
        end

        myfile = fullfile(output, fname);
        fid = fopen(myfile,'w');

        fprintf(fid,'onset\tduration\ttrial_type\tresponse_time\tEndowment\tDecision\tOffer\n');

        %% ---- WRITE EVENTS ----
        for t = 1:length(Block)

            % Cue phase
            if Block(t) == 3
                trial_type = 'cue_social';
            elseif Block(t) == 2
                trial_type = 'cue_nonsocial';
            else
                trial_type = '';
            end

            if ~isempty(trial_type)
                fprintf(fid,'%f\t2\t%s\t%s\t%d\t%s\t%s\n', ...
                    onset(t), trial_type, 'n/a', Endowment(t), 'n/a', 'n/a');

                if Endowment(t) > 20
                    fprintf(fid,'%f\t2\t%s_high\t%s\t%d\t%s\t%s\n', ...
                        onset(t), trial_type, 'n/a', Endowment(t), 'n/a', 'n/a');
                else
                    fprintf(fid,'%f\t2\t%s_low\t%s\t%d\t%s\t%s\n', ...
                        onset(t), trial_type, 'n/a', Endowment(t), 'n/a', 'n/a');
                end
            end

            % Decision phase
            if response(t) == 999
                fprintf(fid,'%f\t3.7669463\tmissed_trial\t%s\t%d\t%s\t%s\n', ...
                    decision_onset(t), 'n/a', Endowment(t), 'n/a', 'n/a');
            else

                if Block(t) == 3 && Endowment(t) > 20
                    trial_type = 'dec_social_high';
                elseif Block(t) == 3 && Endowment(t) < 20
                    trial_type = 'dec_social_low';
                elseif Block(t) == 2 && Endowment(t) > 20
                    trial_type = 'dec_nonsocial_high';
                else
                    trial_type = 'dec_nonsocial_low';
                end

                fprintf(fid,'%f\t%f\t%s\t%f\t%d\t%s\t%d\n', ...
                    decision_onset(t), RT(t), trial_type, RT(t), ...
                    Endowment(t), 'n/a', offer(t));

                if (response(t) == 1 && L_Option(t) > 0) || ...
                   (response(t) == 2 && R_Option(t) > 0)

                    fprintf(fid,'%f\t%f\t%s_accept\t%f\t%d\taccept\t%d\n', ...
                        decision_onset(t), RT(t), trial_type, RT(t), ...
                        Endowment(t), offer(t));
                else
                    fprintf(fid,'%f\t%f\t%s_reject\t%f\t%d\treject\t%d\n', ...
                        decision_onset(t), RT(t), trial_type, RT(t), ...
                        Endowment(t), offer(t));
                end
            end
        end

        fclose(fid);
        fprintf('Created: %s\n', myfile);

    end
end
