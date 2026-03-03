%% This script currently only works from the linux box due to hardcoded paths 

%function convertTrust_BIDS(sublist(s))
%% set up dirs
scriptname = matlab.desktop.editor.getActiveFilename;
[codedir,~,~] = fileparts(scriptname);
% Eliminate cd
cd(codedir);
addpath(codedir)
cd ..
maindir = '/ZPOOL/data/projects/rf1-betrayal';%pwd;

% Partner is Friend=3, Stranger=2, Computer=1
% Reciprocate is Yes=1, No=0
% cLeft is the left option
% cRight is the right option
% high/low value option will randomly flip between left/right
% 
% sublist(s) = sublist(s)ect number
%/data/projects/rf1-sra/stimuli/Scan-Investment_Game/logs
rawdata='/ZPOOL/data/projects/rf1-sra/stimuli/Scan-Investment_Game/logs';

%% ---- LOAD SUBJECT LIST ----
sublist_file = fullfile(maindir, 'code', 'subject_list_n225.txt');

assert(isfile(sublist_file), 'Subject list file not found.')


%% ---- LOAD SUBJECT LIST ----
sublist_raw = readlines(sublist_file);   % read as strings
sublist_raw = strtrim(sublist_raw);      % remove whitespace
sublist_raw = sublist_raw(sublist_raw ~= ""); % remove empty lines

% Remove 'sub-' prefix
sublist_clean = erase(sublist_raw, "sub-");

% Convert to numeric
sublist = str2double(sublist_clean);

% Remove any rows that failed conversion
sublist = sublist(~isnan(sublist));



for s = 1:length(sublist)

try
    
    for r = 0:1
        %r=0;
        fname = fullfile(rawdata,num2str(sublist(s)),sprintf('sub-%03d_task-trust_run-%d_raw.csv',sublist(s),r));
        if exist(fname,'file')
            disp("First checkpoint");
            fid = fopen(fname,'r');
        else
            fprintf('sub-%d -- Investment Game, Run %d: No data found.\n', sublist(s), r+1)
            continue;
        end
        
        %Todo: change from textscan to readtable
        % C = textscan(fid,[repmat('%f',1,14) '%s' repmat('%f',1,9)],'Delimiter',',','HeaderLines',1,'EmptyValue', NaN);
        % fclose(fid);
        % 
        % outcomeonset = C{20}; % should be locked to the presentation of the partner cue (at least 500 ms before choice screen)
        % choiceonset = C{11}; % should be locked to the presentation of the partner cue (at least 500 ms before choice screen)
        % RT = C{18};
        % Partner = C{5};
        % reciprocate = C{6};
        % response = C{15}; % high/low -- build in check below to check recording
        % trust_val = C{13}; % 0-8 (with '999' for no response)
        % cLeft = C{3};
        % cRight = C{4};
        % options = [cLeft cRight];

                % Read the data using readtable, delete textscan later 
        T = readtable(fname, 'Delimiter', ',');
        fclose(fid);

        % Extract the columns you need using column names
        outcomeonset = T.outcome_onset; % should be locked to the presentation of the partner cue (at least 500 ms before choice screen)
        choiceonset = T.onset; % should be locked to the presentation of the partner cue (at least 500 ms before choice screen)
        RT = T.rt;
        Partner = T.Partner;
        reciprocate = T.Reciprocate;
        response = T.highlow; % high/low -- build in check below to check recording
        trust_val = T.resp; % 0-8 (with '999' for no response)
        cLeft = T.cLeft;
        cRight = T.cRight;
        options = [cLeft, cRight];

        
        fname = sprintf('sub-%03d_task-trust_run-%01d_events.tsv',sublist(s),r+1);
        output = fullfile(maindir,'bids',['sub-' num2str(sublist(s))],'func');
        disp(output);
        if ~exist(output,'dir')
            disp("Second checkpoint");
            mkdir(output)
        end
        fid = fopen(fullfile(output,fname),'w');
        fprintf(fid,'onset\tduration\ttrial_type\tresponse_time\ttrust_value\tchoice\tcLow\tcHigh\n');
        for t = 1:length(choiceonset)
            
            % output check
            if trust_val(t) == 999
                if strcmp(response{t},'high') && (max(options(t,:)) ~= trust_val(t))
                    error('response output incorrectly recorded for trial %d', t)
                end
            end
            
            if (Partner(t) == 1)
                trial_type = 'computer';
            elseif (Partner(t) == 2)
                trial_type = 'stranger';
            elseif (Partner(t) == 3)
                trial_type = 'friend';
            end
            
            % "String values containing tabs MUST be escaped using double quotes.
            % Missing and non applicable values MUST be coded as "n/a"."
            % http://bids.neuroimaging.io/bids_spec.pdf
            
            %fprintf(fid,'onset\tduration\ttrial_type\tresponse_time\ttrust_value\tchoice\n');
            if trust_val(t) == 999
                fprintf(fid,'%f\t%f\t%s\t%f\t%s\t%s\t%d\t%d\n',choiceonset(t),3,'missed_trial',3,'n/a','n/a',min(options(t,:)),max(options(t,:)));
            else
                if trust_val(t) == 0
                    fprintf(fid,'%f\t%f\t%s\t%f\t%d\t%s\t%d\t%d\n',choiceonset(t),RT(t),['choice_' trial_type ],RT(t),0,response{t},min(options(t,:)),max(options(t,:))); %should always be 'low'
                else
                    if reciprocate(t) == 1
                        fprintf(fid,'%f\t%f\t%s\t%f\t%d\t%s\t%d\t%d\n',choiceonset(t),RT(t),['choice_' trial_type],RT(t),trust_val(t),response{t},min(options(t,:)),max(options(t,:)));
                        fprintf(fid,'%f\t%f\t%s\t%f\t%d\t%s\t%d\t%d\n',outcomeonset(t),1,['outcome_' trial_type '_recip'],RT(t),trust_val(t),response{t},min(options(t,:)),max(options(t,:)));
                    else
                        fprintf(fid,'%f\t%f\t%s\t%f\t%d\t%s\t%d\t%d\n',choiceonset(t),RT(t),['choice_' trial_type],RT(t),trust_val(t),response{t},min(options(t,:)),max(options(t,:)));
                        fprintf(fid,'%f\t%f\t%s\t%f\t%d\t%s\t%d\t%d\n',outcomeonset(t),1,['outcome_' trial_type '_defect'],RT(t),trust_val(t),response{t},min(options(t,:)),max(options(t,:)));
                    end
                end
            end
            
        end
        fclose(fid);
        % rand_trial = randsample(1:36,1);
        % if trust_val(rand_trial) == 999
        %     %fprintf('sub-%d -- Investment Game, Run %d: On trial %d, Participant did not respond.\n',sublist(s), r+1, rand_trial);
        % else
        %     if reciprocate(rand_trial)
        %         participant = (8 - trust_val(rand_trial)) + ((trust_val(rand_trial) * 3)/2);
        %         friend = (trust_val(rand_trial) * 3)/2;
        %     else
        %         participant = 8 - trust_val(rand_trial);
        %         friend = (trust_val(rand_trial) * 3);
        %     end
        %     if (Partner(rand_trial) == 1)
        %         trial_type = 'Computer';
        %     elseif (Partner(rand_trial) == 2)
        %         trial_type = 'Stranger';
        %     elseif (Partner(rand_trial) == 3)
        %         trial_type = 'Friend';
        %     end
        %     %fprintf('sub-%d -- Investment Game, Run %d: On trial %d, Participant WINS $%.2f and %s WINS $%.2f.\n', sublist(s), r+1, rand_trial, participant, trial_type, friend);
        % end
    end
    
catch ME
    disp(ME.message)
    msg = sprintf('check line %d', ME.stack.line);
    disp(msg);
    keyboard
end
end