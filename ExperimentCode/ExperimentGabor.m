function GaborData = ExperimentGabor(GaborData, varargin)

directory = fullfile(pwd, '..');
settings = LoadSettings(directory);

datadir = fullfile(directory, 'RawData');
if ~exist(datadir, 'dir'), mkdir(datadir); end

subjectID = getSubjectId(datadir, 'gaborV2');
sessionNo = length(dir(fullfile(datadir, [subjectID '*']))) + 1;
subjectID = [subjectID '-Session' num2str(sessionNo)];

%% Environment and PsychToolBox Initialization

if isempty(GaborData.model_observer)
    
    % Define variables that PTB adds to the 'static workspace' here (avoids an
    % error message...)
    global AGL GL GLU ptb_RootPath ptb_ConfigPath;
    
    cd(fullfile(directory, 'Code')) % Set the current directory
    commandwindow; % Moves the cursor to the commandwindow
    
    if settings.useOpenGL, InitializeMatlabOpenGL; end
    
    % Screen set up
    whichScreen = 0;%settings.whichScreen; %allow to choose the display if there's more than one
    

    xc = settings.screenSize(3)/2; %	Gets the middle of the horizontal axis
    yc = settings.screenSize(4)/2; % Gets the middle of the vertical axis
    GaborData.lims = [0+50 0+50 settings.screenSize(3)-50 settings.screenSize(4)-50];
    GaborData.center = [xc yc];
    Screen('Preference', 'SkipSyncTests', settings.ptbSkipSyncTests); % Opens Screen
    
    white = [255 255 255];          % Sets the color to be white
    black = [0 0 0];                % Sets the color to be black
    
    [wPtr, ~] = Screen('OpenWindow', whichScreen, black, [1920 0 1920*2 1080], 32); % Opens window, sets background as black, sets screensize

    
    if ~isempty(settings.gammaTableFile)
        gtdata = load(settings.gammaTableFile);
        Screen('LoadNormalizedGammaTable', wPtr, gtdata.(settings.gammaTable)*[1 1 1]);
    end
    
    
    % Set up keyboard functions
    KbName('UnifyKeyNames');
    goKey = KbName(settings.keyGo);
    exitKey = KbName(settings.keyExit);
    
end

if GaborData.no_staircase
    fileName = fullfile(datadir, [subjectID '-GaborDataNoStaircase.mat']);
    fileNameQuit = fullfile(datadir, [subjectID '-GaborDataNoStaircaseQuit.mat']);
else
    if isequal(GaborData.stair_fn, @Staircase.contrast)
        fileName = fullfile(datadir, [subjectID '-GaborDataContrast.mat']);
        fileNameQuit = fullfile(datadir, [subjectID '-GaborDataContrastQuit.mat']);
    elseif isequal(GaborData.stair_fn, @Staircase.ratio)
        fileName = fullfile(datadir, [subjectID '-GaborDataRatio.mat']);
        fileNameQuit = fullfile(datadir, [subjectID '-GaborDataRatioQuit.mat']);
    elseif isequal(GaborData.stair_fn, @Staircase.noise)
        fileName = fullfile(datadir, [subjectID '-GaborDataNoise.mat']);
        fileNameQuit = fullfile(datadir, [subjectID '-GaborDataNoiseQuit.mat']);
    else
        warning('No staircase-specific suffix on save name!');
        fileName = fullfile(datadir, [subjectID '.mat']);
        fileNameQuit = fullfile(datadir, [subjectID 'Quit.mat']);
    end
end

if exist(fileName, 'file')
    Screen('CloseAll');
    error('Data for %s already exists', fileName);
end

[tracker_info] = Eyelink.Initialize_params(whichScreen,wPtr,'eyelink_use', GaborData.eyelink_use);
tracker_info = Eyelink.setup(tracker_info);%,fileName_edf,edfdir);

if tracker_info.eyelink_use
    HideCursor(whichScreen);
end

if isempty(GaborData.model_observer)
    % Create 2 textures to display templates.
    right_template = squeeze(bpg.genImages(1, GaborData.stim_size, GaborData.stim_sp_freq_cpp, GaborData.stim_std_sp_freq_cpp, GaborData.right_category, 2)) * 64.0 + 127.0;
    left_template = squeeze(bpg.genImages(1, GaborData.stim_size, GaborData.stim_sp_freq_cpp, GaborData.stim_std_sp_freq_cpp, GaborData.left_category, 2)) * 64.0 + 127.0;
    right_tex = Screen('MakeTexture', wPtr, right_template);
    left_tex = Screen('MakeTexture', wPtr, left_template);
    [h, w] = size(right_template);
    
    
    %% Begin experiment
    if tracker_info.eyelink_use
        tracker_info = Eyelink.calibrate(tracker_info);
    end
    
    % Instruction Screen
    textbox = ptbCenteredRect([xc yc], settings.screenSize(3:4)/3);
    Screen('TextSize', wPtr, 20); % Set text size to 30
    Screen('FillRect', wPtr, 127);
    
    
    DrawFormattedText(wPtr, ...`
        ['You will see a series of images flashing very quickly. ' ...`
        'You are required to keep your eyes on the white cross. While fixating on a pattern stimulus, you will be shown two equally-distant patterns in the periphery. ' ...`
        'You will have to choose to saccade to one of the two patterns in the periphery; and make a quick eyemovement to the chosen pattern once the fixated pattern disappears. This would help you decide the correct overall pattern of the frames. ' ...`
        'At the end of the trial, you will have to choose the pattern that is more consistent with all of the preceding frames (as shown below). ' ...`
        sprintf('Select the pattern positioned to the left or right by pressing %s or %s respectively. ', settings.keyLeftName, settings.keyRightName) ...`
        'Ask the researcher if you need further clarification. ' ...`
        sprintf('Press %s to begin.', settings.keyGoName)], ...`
        'centerblock', 'center', white, 100, 0, 0, 1.5, 0, textbox);
    
    
    
    
    templatey = (textbox(4) + settings.screenSize(4)) /2;
    Screen('DrawTexture', wPtr, left_tex, [], ptbCenteredRect([xc-w templatey], [w h]));
    Screen('DrawTexture', wPtr, right_tex, [], ptbCenteredRect([xc+w templatey], [w h]));
    Screen('Flip', wPtr); % Function to flip to the next screen image
    if ptbWaitKey([goKey exitKey]) == exitKey
        Eyelink.finish(tracker_info);
        Screen('CloseAll');
        return;
    end
    Screen('Flip', wPtr); % Function to flip to the next screen image
end

    function earlyQuit
        if GaborData.current_trial > 5
            save(fileNameQuit, 'GaborData');
        end
        ShowCursor();
        Eyelink.finish(tracker_info);
        Screen('CloseAll');
    end

% Begin Preliminary Trials
seen_block_notification = false;
trial = 1;
block_trial = 1;
block = 1;

% Using 'while' rather than 'for' since invalid trials (broke fixation or
% didn't respond in time) don't increment 'trial'.
try
    while trial <= GaborData.trials_per_block * GaborData.blocks
        
        %% Bookkeeping to set up trial
        
        GaborData.current_trial = trial;
        %just for preliminary stimulus generation
        
        % Reset params at the start of each block
        if mod(trial, GaborData.trials_per_block) == 1
            % Display a message if this is the beginning of the second or
            % higher block.
            if trial ~= 1 && isempty(GaborData.model_observer) && ~seen_block_notification
                seen_block_notification = true;
                if isempty(GaborData.model_observer)
                    sounds(-1, 1.5);
                    DrawFormattedText(wPtr, ...
                        [sprintf('You have completed a %d blocks. ', block) ...
                        'You may take a break if you want! ' ...
                        sprintf('Press %s whenever you are ready again.', settings.keyGoName) ...
                        '\nThe images to be discriminated are displayed again below.'], ...
                        'centerblock', 'center', white, 60, 0, 0, 1.5, 0, textbox);
                    Screen('DrawTexture', wPtr, left_tex, [], ptbCenteredRect([xc-w templatey], [w h]));
                    Screen('DrawTexture', wPtr, right_tex, [], ptbCenteredRect([xc+w templatey], [w h]));
                    Screen('Flip', wPtr);
                    block=block+1;
                    save(fileName, 'GaborData');
                    if ptbWaitKey([goKey exitKey]) == exitKey
                        earlyQuit;
                        return;
                    end
                    
                    Screen('Flip', wPtr);
                     if tracker_info.eyelink_use
        tracker_info = Eyelink.calibrate(tracker_info);
    end
                end
            end
            
            % Start of a block - set params to initial values.
            GaborData.streak(trial) = 0;
            GaborData.reversal_counter(trial) = 0;
            GaborData.contrast(trial) = GaborData.contrast(1);
            GaborData.ratio(trial) = GaborData.ratio(1);
            GaborData.noise(trial) = GaborData.noise(1);
            GaborData.step_size(trial) = GaborData.step_size(1);
            if isfield(GaborData, 'iid')
                GaborData.iid(trial) = GaborData.iid(1);
            end
            block_trial = 1;
        else
            seen_block_notification = false;
            
            GaborData.contrast(trial) = GaborData.contrast(trial-1);
            GaborData.ratio(trial) = GaborData.ratio(trial-1);
            GaborData.noise(trial) = GaborData.noise(trial-1);
            GaborData.step_size(trial) = GaborData.step_size(trial-1);
            
            % Count correct streak (with respect to the ideal observer's
            % answer, not the underlying distribution)
            if GaborData.ideal_answer(trial-1) == GaborData.choice(trial-1)
                GaborData.streak(trial) = GaborData.streak(trial-1) + 1;
            else
                GaborData.streak(trial) = 0;
            end
            
            % Count reversals
            if block_trial > 2 && sign(GaborData.streak(trial-1)) ~= sign(GaborData.streak(trial))
                GaborData.reversal_counter(trial) = GaborData.reversal_counter(trial-1) + 1;
            else
                GaborData.reversal_counter(trial) = GaborData.reversal_counter(trial-1);
            end
            
            % Apply the staircase
            if ~GaborData.no_staircase
                GaborData = GaborData.stair_fn(GaborData);
            else
                GaborData.contrast(trial) = GaborData.contrast(trial-1);
                GaborData.noise(trial) = GaborData.noise(trial-1);
                GaborData.ratio(trial) = GaborData.ratio(trial-1);
            end
            
        end
        % Run this trial
        
        % Generate stimulus for this trial.
        [image_array, frame_categories, checksum] = GaborStimulus(GaborData, trial);
        %         GaborData.frame_categories(trial, :,:) = frame_categories;
        GaborData.frame_categories{trial} = frame_categories;
        GaborData.checksum(trial) = checksum;
        image_array_for_sig = [];
        image_array_new = [];
        if isempty(GaborData.model_observer)
            % Pass in the contrast level, all of the images, screen being
            % used, subject ID, the struct with all of the data, and the
            % fact it's the person or computer running the experiment
            image_array_new = reshape(image_array,GaborData.num_peri(trial),GaborData.number_of_images+1,GaborData.stim_size,GaborData.stim_size);
            for pr=2:GaborData.num_peri(trial)
                image_array_new(pr,1,:,:) = image_array_new(1,1,:,:);
            end
            %to compute ideal answer we store images on screen
            image_array_for_sig(1,:,:) = squeeze(image_array_new(1,1,:,:));
            tmp_indx = 2;
            for fr=2:GaborData.number_of_images+1
                for pr=1:GaborData.num_peri(trial)
                    image_array_for_sig(tmp_indx,:,:) = squeeze(image_array_new(pr,fr,:,:));
                    tmp_indx = tmp_indx + 1;
                end
            end
            [I, locations, tracker_info,broke_condition, quit] = trialStimuliGabor(GaborData, unit8(image_array_new), wPtr, tracker_info, settings);
            %GaborData.image_array_new{trial}=image_array_new;
            GaborData.stim_locations{trial} = locations;
            GaborData.ideal_frame_signals{trial} = ...
                bpg.getSignal(double(image_array_for_sig) - 127, GaborData.left_category, max(GaborData.noise(trial), .04)) - ...
                bpg.getSignal(double(image_array_for_sig) - 127, GaborData.right_category, max(GaborData.noise(trial), .04));
            
            GaborData.notchosen_ideal_frame_signals{trial} = ...
                bpg.getSignal(double(I.image_array_notchosen) - 127, GaborData.left_category, max(GaborData.noise(trial), .04)) - ...
                bpg.getSignal(double(I.image_array_notchosen) - 127, GaborData.right_category, max(GaborData.noise(trial), .04));
            
            GaborData.chosen_ideal_frame_signals{trial} = ...
                bpg.getSignal(double(I.image_array_chosen) - 127, GaborData.left_category, max(GaborData.noise(trial), .04)) - ...
                bpg.getSignal(double(I.image_array_chosen) - 127, GaborData.right_category, max(GaborData.noise(trial), .04));
            
            GaborData.ideal_answer(trial) = 1 * (sum(GaborData.ideal_frame_signals{trial}) > 0);
            %GaborData.image_array_chosen{trial} = I.image_array_chosen;
            %GaborData.image_array_notchosen{trial} = I.image_array_notchosen;
            
            GaborData.image_array_chosen_index{trial} = I.image_indx_chosen;
            GaborData.image_array_notchosen_index{trial} = I.image_indx_notchosen;
            
            GaborData.total_number_of_images(trial) = GaborData.number_of_images * GaborData.num_peri(trial) + 1;
            %GaborData.true_ratio(trial) = (sum(GaborData.frame_categories{trial}(:,2:end),'all') + GaborData.frame_categories{trial}(1,1))/45/(GaborData.total_number_of_images(trial));
            
            true_ratio_right = (length(find(GaborData.frame_categories{trial}(:,2:end)==45)) + logical(GaborData.frame_categories{trial}(1,1)==45))/(GaborData.total_number_of_images(trial));
            true_ratio_left = (length(find(GaborData.frame_categories{trial}(:,2:end)==-45)) + logical(GaborData.frame_categories{trial}(1,1)==-45))*(-1)/(GaborData.total_number_of_images(trial));
            
            if abs(true_ratio_right)>abs(true_ratio_left)
                GaborData.true_ratio(trial)=true_ratio_right;
            else
                GaborData.true_ratio(trial)=true_ratio_left;
            end
            
            % Add a 'signed noise' and 'signed contrast' field.
            trial_sign = sign(GaborData.true_ratio(trial) - 0.5);
            GaborData.sign_noise(trial) = trial_sign .* GaborData.noise(trial);
            GaborData.sign_contrast(trial) = trial_sign .* GaborData.contrast(trial);
            
            
            if quit, earlyQuit; return; end
            
            if broke_condition || isnan(I.choice)
                
                Screen('FillRect', wPtr, 127);
                Screen('Flip', wPtr);
                %sounds(2, 0.2);
                WaitSecs(1);
                continue;
            end
            
            GaborData.choice(trial) = I.choice;
            GaborData.reaction_time(trial) = I.reaction;
            GaborData.eye_tracker_points{trial} = tracker_info;
            GaborData.gaze_choice{trial} = I.gazechoice;
        elseif strcmpi(GaborData.model_observer, 'ideal')
            GaborData.choice(trial) = GaborData.ideal_answer(trial);
        elseif strcmpi(GaborData.model_observer, 'oracle')
            GaborData.choice(trial) = GaborData.correct_answer(trial);
        elseif strcmpi(GaborData.model_observer, 'bernoulli')
            decision_var = dot(GaborData.ideal_frame_signals(trial, :), GaborData.model_pk);
            bernoulli_p = sigmoid(decision_var / GaborData.sigmoid_slope);
            % Choose sign of decision_var with probability related to
            % magnitude of decision_var.
            GaborData.choice(trial) = 1 * (rand < bernoulli_p);
        end
        
        %% Accuracy & Feedback
        GaborData.accuracy(trial) = GaborData.choice(trial) == GaborData.correct_answer(trial);
        if isempty(GaborData.model_observer)
      stimulus_bbox = ptbCenteredRect([settings.screenSize(3)/2 settings.screenSize(4)/2], [100 100]);

            if GaborData.accuracy(trial)
             %   sounds(1, 0.2);

            Screen('FillOval',wPtr,[0 0 0],stimulus_bbox);

            else
             %   sounds(0, 0.2);
            Screen('FillOval',wPtr,[255 255 255],stimulus_bbox);

            end
            Screen('Flip', wPtr);
            WaitSecs(0.5); % Pause for 500 ms after feedback before next trial
        end
        
        trial = trial + 1;
        block_trial = block_trial + 1;
    end
    
catch ERR
    earlyQuit;
    Screen('CloseAll');
    rethrow(ERR);
end

%% Save final data to folder
Screen('CloseAll');
Eyelink.finish(tracker_info);
save(fileName, 'GaborData');
end
