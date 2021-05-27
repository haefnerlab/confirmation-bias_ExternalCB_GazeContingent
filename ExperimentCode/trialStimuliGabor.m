function [image_properties, location, tracker_info, broke_condition, quit] = trialStimuliGabor(Data, image_array, wPtr, tracker_info, settings)
% trialStimuli displays the animation of several gabor patches in quick
% succession to the subject, or runs through a single trial of the
% experiment.
if tracker_info.eyelink_use
    tracker_info = Eyelink.startTrialPrepare(tracker_info);
end
image_properties = struct();
quit = false;
broke_condition = false;


%% Make sure to have left/right patch to match the orientations used

% Create images to be displayed as left or right options
left_patch = squeeze(bpg.genImages(1, Data.stim_size, Data.stim_sp_freq_cpp, Data.stim_std_sp_freq_cpp, Data.left_category, 2)) * 64.0 + 127.0;
right_patch = squeeze(bpg.genImages(1, Data.stim_size, Data.stim_sp_freq_cpp, Data.stim_std_sp_freq_cpp, Data.right_category, 2)) * 64.0 + 127.0;

xc = settings.screenSize(3)/2; % Get the side of the horizontal axis
yc = settings.screenSize(4)/2; % Get the side of the vertical axis

% black = [0 0 0];
gray = [127 127 127];

% Set up variables for keyboard functions
KbName('UnifyKeyNames');
exitKey = KbName(settings.keyExit);
leftKey = KbName(settings.keyLeft);
rightKey = KbName(settings.keyRight);

total_frames = Data.number_of_images;
%Getting all the possible location of stimuli for all the frames.
location = getLocation(Data.lims, Data.center, Data.number_of_images+1, Data.num_peri(Data.current_trial)); % number_of_images+1 to be able to have locations for periphery of the lat foveated frame
image_texture = zeros(Data.num_peri(Data.current_trial), total_frames+1);
%Getting the image texture for all the frames
for i = 1:Data.number_of_images+1
    for j=1:Data.num_peri(Data.current_trial)
        image_texture(j,i) = Screen('MakeTexture', wPtr, squeeze(image_array(j,i, :, :)));
    end
end
[~, ~,h, w] = size(image_array);
image_properties.image_array_chosen = zeros(total_frames,h,w);
image_properties.image_array_notchosen = zeros(((Data.num_peri(Data.current_trial)-1) * (total_frames-1) + Data.num_peri(Data.current_trial)),h,w); % beacuse in the last foveated frame two images in the periphery were shown but not chosen
show_left_patch = Screen('MakeTexture', wPtr, left_patch);
show_right_patch = Screen('MakeTexture', wPtr, right_patch);


Screen('FillRect', wPtr, gray);  % Make the background gray
Screen('Flip', wPtr);


% drawTrialNo();

% pos refers to the choice made in the current frame
% selected im refers to selected image in the previous frame.
pos = 1;
posp = [];
selected_im = 1;
% Preparing for eye-tracking protocol when eyelink_use is true.
if tracker_info.eyelink_use
    tracker_info = Eyelink.startTrial(tracker_info);
    Eyelink.clearBuffer(tracker_info.drained);
end
image_properties.image_indx_chosen = false(total_frames+1,Data.num_peri(Data.current_trial));
image_properties.image_indx_notchosen = false(total_frames+1,Data.num_peri(Data.current_trial));

temp_length=0;
n_ind = 1;
image_properties.gazechoice(1,:) = location{1};
image_properties.image_indx_chosen(1,:) = true(1,Data.num_peri(Data.current_trial));
image_properties.image_indx_notchosen(1,:) = false(1,Data.num_peri(Data.current_trial));
image_properties.image_array_chosen(1,:,:) = image_array(1,1,:,:);
disp('Stim Starts');
for i = 1:total_frames
    % fixation reference
    
    if i==1
        EyeTracker.drawFixationSymbol(tracker_info,location{i}(pos,1),location{i}(pos,2), wPtr);
        [start_time,~] = Screen('Flip', wPtr);
        EyeTracker.drawFixationSymbol(tracker_info,location{i}(pos,1),location{i}(pos,2), wPtr);
        [first_fixation,~] = Screen('Flip', wPtr, start_time + Data.settle_time);
        % Tracking the eye position
        [tracker_info, broke_fixation, ~, temp_length] = Eyelink.track_eye(tracker_info,temp_length,location{i}(pos,:), Data.first_fixation_duration, 1,'is_fixating');
        % If the fixation broke, the trial aborts. Otherwise, continues.
        if broke_fixation==1
            broke_condition = 1;
            image_properties.choice=nan;
            return;
        end
        % Only for the first frame, allowing more time (first_fixation_duration) to stablize the
        % fixation on the main stimulus.
        time_flip = first_fixation + Data.first_fixation_duration;
        
        stimulus_box = ptbCenteredRect([location{i}(pos,1),location{i}(pos,2)], [w h]);
        Screen('DrawTexture', wPtr, image_texture(selected_im,i), [], stimulus_box); %Fill the buffer with the first texture
        EyeTracker.drawFixationSymbol(tracker_info,location{i}(pos,1),location{i}(pos,2), wPtr);
        % Presenting the first main stimulus.
        [stim_onset,~] = Screen('Flip', wPtr, time_flip);
        
        [tracker_info, broke_fixation, ~, temp_length ] = Eyelink.track_eye(tracker_info, temp_length, location{i}(pos,:), Data.pure_frame_duration, 1, 'is_fixating');
        if broke_fixation==1
            broke_condition = 1;
            image_properties.choice=nan;
            return;
        end
        
    else
        time_flip = 0;
        stimulus_box = ptbCenteredRect([location{i}(pos,1),location{i}(pos,2)], [w h]);
        Screen('DrawTexture', wPtr, image_texture(selected_im,i), [], stimulus_box); %Fill the buffer with the first texture
        EyeTracker.drawFixationSymbol(tracker_info,location{i}(pos,1),location{i}(pos,2), wPtr);
        % Presenting the first main stimulus.
        [stim_onset,~] = Screen('Flip', wPtr, time_flip);
        
        [tracker_info, broke_fixation, ~, temp_length ] = Eyelink.track_eye(tracker_info, temp_length, location{i}(pos,:), Data.pure_frame_duration, 1, 'is_fixating');
        if broke_fixation==1
            broke_condition = 1;
            image_properties.choice=nan;
            return;
        end
    end
    for pr=1:Data.num_peri(Data.current_trial)
        posp(pr) = (pos-1) * Data.num_peri(Data.current_trial) + pr;
    end
    
    % for all the frames except for the last one, defining references for two possible
    % locations to land saccade (pos1 and pos2).
    
    % Presenting main stimulus + two randomly selected periphery
    % stimuli.
    Screen('DrawTexture', wPtr, image_texture(selected_im,i), [], stimulus_box);
    for pr=1:Data.num_peri(Data.current_trial)
        stimulus_box_peri(pr,:) = ptbCenteredRect([location{i+1}(posp(pr),1),location{i+1}(posp(pr),2)], [w h]);
        Screen('DrawTexture', wPtr, image_texture(pr,i+1), [], stimulus_box_peri(pr,:));
        %EyeTracker.drawFixationSymbol(tracker_info,location{i+1}(posp(pr),1),location{i+1}(posp(pr),2), wPtr);
    end
    EyeTracker.drawFixationSymbol(tracker_info,location{i}(pos,1),location{i}(pos,2), wPtr);
    [peri_onset, ~] = Screen('Flip', wPtr, stim_onset + Data.pure_frame_duration);
    
    [tracker_info,broke_fixation, ~, temp_length] = Eyelink.track_eye(tracker_info, temp_length, location{i}(pos,:), Data.periphery_duration, 1, 'is_fixating');
    if broke_fixation==1
        broke_condition = 1;
        image_properties.choice=nan;
        return;
    end
    if i<(total_frames)
        % Removing the foveated stimulus and keeping onlythe periphery
        % stimuli to cue saccade initiation
        for pr=1:Data.num_peri(Data.current_trial)
            Screen('DrawTexture', wPtr, image_texture(pr,i+1), [], stimulus_box_peri(pr,:));
            EyeTracker.drawFixationSymbol(tracker_info,location{i+1}(posp(pr),1),location{i+1}(posp(pr),2), wPtr);
        end
            EyeTracker.drawFixationSymbol(tracker_info,location{i}(pos,1),location{i}(pos,2), wPtr);

        Screen('Flip', wPtr, peri_onset + Data.periphery_duration);
        
        % The track_eye function also returns the location where saccade
        % landed (saccade_loc).
        temp_loc = [];
        for pr=1:Data.num_peri(Data.current_trial)
            temp_loc = [temp_loc; location{i+1}(posp(pr),:)];
        end
        [tracker_info, not_saccaded, saccade_loc, temp_length] = Eyelink.track_eye(tracker_info, temp_length, temp_loc, Data.saccade_limit,1,'saccading');
        if not_saccaded==1
            broke_condition = 1;
            image_properties.choice=nan;
            return;
        else
            % If the saccade was made, figuring out where the saccade
            % landed by comparing with the two possible location using the
            % location reference defined earlier.
            for pr=1:size(temp_loc,1)
                if saccade_loc == temp_loc(pr,:)
                    selected_im = pr;
                    pos = posp(pr);
                    rem = mod(pr,Data.num_peri(Data.current_trial));
                end
            end
            if rem==0
                rem = Data.num_peri(Data.current_trial);
            end
            not_selected_im = setdiff([1:Data.num_peri(Data.current_trial)],rem);
            temp_yes = false(1,Data.num_peri(Data.current_trial));
            temp_yes(rem) = true;
            image_properties.image_indx_chosen(i+1,:) = temp_yes; 
            temp_not = false(1,Data.num_peri(Data.current_trial));
            temp_not(not_selected_im) = true;
            image_properties.image_indx_notchosen(i+1,:) = temp_not;
            image_properties.gazechoice(i+1,:) = temp_loc(selected_im,:,:);
            image_properties.image_array_chosen(i+1,:,:) = image_array(rem,i+1,:,:);
            if Data.num_peri(Data.current_trial)>1
                image_properties.image_array_notchosen(n_ind:n_ind+length(not_selected_im)-1,:,:) = image_array(not_selected_im,i+1,:,:);
            end
            n_ind = n_ind+length(not_selected_im);
        end
    end
end

for pr=1:Data.num_peri(Data.current_trial)
    image_properties.image_array_notchosen(n_ind,:,:) = image_array(pr,end,:,:);
    n_ind = n_ind +1;
end
image_properties.image_indx_notchosen(end,:) = true(1,Data.num_peri(Data.current_trial));
image_properties.image_indx_chosen(end,:) = false(1,Data.num_peri(Data.current_trial));

Screen('FillRect', wPtr, gray);
[~, endTime] = Screen('Flip', wPtr);
Screen('Flip', wPtr, endTime + Data.go_cue_time);

Screen('DrawTexture', wPtr, show_left_patch, [], ptbCenteredRect([xc-w yc], [w h]));   % xc, yc indicates the coordinates of the middle of the screen
Screen('DrawTexture', wPtr, show_right_patch, [], ptbCenteredRect([xc+w yc], [w h]));
%Screen('DrawText', wPtr, sprintf('Current Trial - #%d', Data.current_trial), xc-600, yc+250, 0);   % Unobtrusive output to screen of the current trial number
Screen('Flip', wPtr);

[key, rt, timeout] = ptbWaitKey([leftKey, rightKey, exitKey], 1);

% Close textures to avoid memory problems.
for i = 1:total_frames
    Screen('Close', image_texture(:,i));
end
Screen('Close', show_left_patch);
Screen('Close', show_right_patch);


if key == exitKey
    quit = true;
    image_properties.choice = nan;
end

if timeout
    image_properties.choice = nan;
else
    image_properties.reaction = rt * 1000;
    if key == leftKey
        image_properties.choice = 1;
    elseif key == rightKey
        image_properties.choice = 0;
    end
end

    function drawTrialNo()
        Screen('DrawText', wPtr, sprintf('Current Trial - #%d', Data.current_trial), xc-900, yc+550, 0);
    end
end