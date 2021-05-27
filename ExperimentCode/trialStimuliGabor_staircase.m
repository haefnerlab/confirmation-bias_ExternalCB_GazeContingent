function [image_properties, tracker_info, broke_condition, quit] = trialStimuliGabor_staircase(Data, image_array, wPtr, tracker_info, settings)
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

xc = settings.screenSize(3)/2; % Get the side of the horizontal axis
yc = settings.screenSize(4)/2; % Get the side of the vertical axis

white = [255 255 255];
black = [0 0 0];
gray = [127 127 127];
location = getLocation(Data.lims,Data.center,Data.number_of_images);

% Set up variables for keyboard functions
KbName('UnifyKeyNames');
exitKey = KbName(settings.keyExit);
leftKey = KbName(settings.keyLeft);
rightKey = KbName(settings.keyRight);

image_texture(1) = Screen('MakeTexture', wPtr, squeeze(image_array(1, :, :)));
image_texture(2) = Screen('MakeTexture', wPtr, squeeze(image_array(2, :, :)));

[~,h, w] = size(image_array);

Screen('FillRect', wPtr, gray);  % Make the background gray
Screen('Flip', wPtr);

% drawTrialNo();

% Preparing for eye-tracking protocol when eyelink_use is true. 
if tracker_info.eyelink_use
    tracker_info = Eyelink.startTrial(tracker_info);
    Eyelink.clearBuffer(tracker_info.drained);
end
temp_length=0;


disp('Stim Starts');
% Presenting fixation symbol for settle time and track eye to make sure the
% gaze is in the location boundary. 
EyeTracker.drawFixationSymbol(tracker_info,location{1}(1),location{1}(2), wPtr);
[start_time,~] = Screen('Flip', wPtr);
EyeTracker.drawFixationSymbol(tracker_info,location{1}(1),location{1}(2), wPtr);
[first_fixation,~] = Screen('Flip', wPtr, start_time + Data.settle_time);
[tracker_info,broke_fixation, ~, temp_length] = Eyelink.track_eye(tracker_info, temp_length, location{1}, Data.first_fixation_duration, 1,'is_fixating');
% If fixation broke, the trial aborts. Otherwise, continues.
if broke_fixation==1
    broke_condition = 1;
    image_properties.choice=nan;
    return;
end
% Presenting the periphery stimulus.
stimulus_box1 = ptbCenteredRect([location{2}(1,1),location{2}(1,2)], [w h]);
stimulus_box2 = ptbCenteredRect([location{2}(2,1),location{2}(2,2)], [w h]);
Screen('DrawTexture', wPtr, image_texture(1), [], stimulus_box1);
Screen('DrawTexture', wPtr, image_texture(2), [], stimulus_box2);
EyeTracker.drawFixationSymbol(tracker_info,location{1}(1,1),location{1}(1,2), wPtr);
% EyeTracker.drawFixationSymbol(tracker_info,location{2}(1,1),location{2}(1,2), wPtr);
% EyeTracker.drawFixationSymbol(tracker_info,location{2}(2,1),location{2}(2,2), wPtr);
[peri_onset, ~] = Screen('Flip', wPtr, first_fixation + Data.first_fixation_duration);
[tracker_info,broke_fixation, ~, temp_length] = Eyelink.track_eye(tracker_info, temp_length, location{1}, Data.periphery_duration, 1, 'is_fixating');
if broke_fixation==1
    broke_condition = 1;
    image_properties.choice=nan;
    return;
end

%Screen('DrawText', wPtr, sprintf('Current Trial - #%d', Data.current_trial), xc-600, yc+250, 0);   % Unobtrusive output to screen of the current trial number


textbox = ptbCenteredRect([xc yc], settings.screenSize(3:4)/3);
Screen('TextSize', wPtr, 20); % Set text size to 30
Screen('FillRect', wPtr, 127);
DrawFormattedText(wPtr, ...
    [sprintf('Select SAME or DIFFERENT by pressing %s or %s respectively. ', settings.keyLeftName, settings.keyRightName)], ...
    'centerblock', 'center', white, 100, 0, 0, 1.5, 0, textbox);
[~, endTime] = Screen('Flip', wPtr, peri_onset + Data.periphery_duration);

Screen('Flip', wPtr, endTime + Data.go_cue_time);

[key, rt, timeout] = ptbWaitKey([leftKey, rightKey, exitKey], 1);

% Close textures to avoid memory problems.

Screen('Close', image_texture(1));
Screen('Close', image_texture(2));


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