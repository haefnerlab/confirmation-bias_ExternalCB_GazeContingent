function GaborData = newGaborData(varargin)

if ~isempty(varargin) && isstruct(varargin{1})
    % Get first input as a 'template' struct.
    template = varargin{1};
    varargin = varargin(2:end);
else
    template = struct();
end

    function value = get_arg(name, default)
        % Helper function to get named arguments with a default
        idx = strcmpi(name, varargin);
        if any(idx)
            val_idx = find(idx)+1;
            value = varargin{val_idx};
            varargin(find(idx):val_idx) = [];
        elseif isfield(template, name)
            if any(strcmpi(name, {'stair_bounds', 'step_size', 'min_step_size'})) && ~isequal(GaborData.stair_fn, template.stair_fn)
                warning('Not copying field %s from template since template''s stair_fn is %s', name, func2str(template.stair_fn));
                value = default;
                return;
            end
            if isnumeric(template.(name))
                expected_size = length(default);
                value = template.(name)(1:expected_size);
            else
                value = template.(name);
            end
        else
            value = default;
        end
    end

%% User-settable params
GaborData.trials_per_block = get_arg('trials_per_block', 50);
GaborData.blocks = get_arg('blocks', 4);
GaborData.stair_fn = get_arg('stair_fn', @Staircase.noise);
GaborData.reversals_per_epoch = get_arg('reversals_per_epoch', 6);
GaborData.no_staircase = get_arg('no_staircase',1);
total_trials = GaborData.trials_per_block * GaborData.blocks;

% Initial values of staircase-able parameters
GaborData.contrast = zeros(1, total_trials);
GaborData.contrast(1) = get_arg('contrast', 60);
GaborData.ratio = zeros(1, total_trials);
GaborData.ratio(1) = get_arg('ratio', .7);
GaborData.noise = zeros(1, total_trials);
GaborData.noise(1) = get_arg('noise', 0.8); % kappa of bpg orientation band
GaborData.step_size = zeros(1, total_trials);

% Staircase bounds and step size, with defaults set depending on stair_fn
GaborData.model_observer = get_arg('model_observer', '');
if isequal(GaborData.stair_fn, @Staircase.contrast)
    GaborData.stair_bounds = get_arg('stair_bounds', [0 64]);
    GaborData.step_size(1) = get_arg('step_size', 2); % multiplicative (in the "easier" direction)
    GaborData.min_step_size = get_arg('min_step_size', 1+(GaborData.step_size(1) - 1)/4); % Default to two 'halvings' of the step size
    GaborData.test_threshold = get_arg('test_threshold', 0);
    %     GaborData.test_ratio = get_arg('test_ratio', 0.9);
    GaborData.test_ratio = get_arg('test_ratio', GaborData.ratio(1));
elseif isequal(GaborData.stair_fn, @Staircase.ratio)
    GaborData.stair_bounds = get_arg('stair_bounds', [0.5 1.0]);
    GaborData.step_size(1) = get_arg('step_size', .1); % additive (in the "easier" direction)
    GaborData.min_step_size = get_arg('min_step_size', GaborData.step_size(1)/4); % Default to two 'halvings' of the step size
elseif isequal(GaborData.stair_fn, @Staircase.noise)
    % Note: noise is treated as a special case, where we use a discrete set
    % of values for kappa, and the staircase simply increments/decrements
    % the index.
    GaborData.kappa_set = get_arg('kappa_set', linspace(0, 0.8, 21)); % Higher indices correspond to easier values of kappa.
    GaborData.stair_bounds = get_arg('stair_bounds', [1 length(GaborData.kappa_set)]); % Not actually used; bounds implied by length of array. See Staircase.noise
    GaborData.step_size(1) = get_arg('step_size', 4); % additive (in the "easier" direction)
    GaborData.min_step_size = get_arg('min_step_size', 1); % Cannot step fewer than 1 indices in an array.
    GaborData.test_threshold = get_arg('test_threshold', 0.4);
    GaborData.test_ratio = get_arg('test_ratio', 0.8);
end
GaborData.no_staircase = 1;
GaborData.eyelink_use = get_arg('eyelink_use', 0);
% Other misc. user-definable parameters relating to stimulus/rig.
GaborData.flag_use_old_stimulus_code = false;  % Henceforth all stimuli are generated using 'correct' code.
GaborData.number_of_images = get_arg('number_of_images', 10);
GaborData.stimulus_fps = get_arg('stimulus_fps', 12);  % frame rate of stimuli
GaborData.blank_frames = get_arg('blank_frames', 5);  % number of blank screen frames per stimulus frame
GaborData.first_fixation_duration = get_arg('first_fixation_duration', 0.2);
GaborData.blank_duration = get_arg('blank_duration',.5);
GaborData.choice_duration = get_arg('choice_duration',2);
%GaborData.choice_timelimit = get_arg('choice_timelimit',1.5);
GaborData.pure_frame_duration = get_arg('pure_frame_duration',0.250);
GaborData.periphery_duration = get_arg('periphery_duration',0.250);
GaborData.saccade_limit = get_arg('saccade_limit',2);
GaborData.settle_time = get_arg('settle_time',.3);

% GaborData.withoutfixation_duration = get_arg('withoutfixation_duration',1);
% GaborData.cue_duration = get_arg('cue_duration', 2);  % Fixed duration, seconds to display cue after getting fixation.
GaborData.annulus = get_arg('annulus', 25); % Size, in pixels, of hole in center of stimulus
GaborData.left_category = get_arg('left_category', +45);
GaborData.right_category = get_arg('right_category', -45);
% GaborData.sanity = get_arg('sanity', 0);
GaborData.go_cue_time = get_arg('go_cue_time', 1);  % Time between final stimulus/mask frame and the targets appearing.

%hexagonal grids
GaborData.lims = [50 50 1870 1030];
GaborData.center = [960 540];
GaborData.randomize_peri = 1;
if GaborData.randomize_peri
    GaborData.num_peri = randi(3,1,total_trials);
else
    GaborData.num_peri = 2 * ones(1,total_trials);
end
% BPG Stimulus parameters
GaborData.stim_size = get_arg('stim_size', 120);  % Width of the stimulus in pixels.
GaborData.stim_sp_freq_cpp = get_arg('stim_sp_freq_cpp', 0.1194);  % Mean spatial frequency of images in cycles per pixel.
GaborData.stim_std_sp_freq_cpp = get_arg('stim_std_sp_freq_cpp', .0597);  % Std deviation of spatial frequency in cycles per pixel.

% Preallocate fields that will be populated with data by running the
% experiment.
GaborData.iid = true(1, total_trials);
GaborData.streak = zeros(1, total_trials);
GaborData.reversal_counter = zeros(1, total_trials);
GaborData.ideal_answer = zeros(1, total_trials);
GaborData.reaction_time = zeros(1, total_trials);
GaborData.choice = zeros(1, total_trials);
GaborData.accuracy = zeros(1, total_trials);

% Note that 'seed' and 'correct_answer' must be preset due to esoteric
% properties of random number generators. GaborStimulus will read out these
% preset values.
GaborData.seed = randi(1000000000, 1, total_trials);
GaborData.checksum = zeros(1, total_trials);  % For sanity-checks on seeds
GaborData.correct_answer = rand(1, total_trials) < .5;



GaborData.current_trial = 0;

GaborData.eye_tracker_points = {};

if isequal(GaborData.stair_fn, @Staircase.ratio)
    if GaborData.step_size(1) < 0
        warning('Changing sign of ratio step_size from %d to %d', GaborData.step_size(1), -GaborData.step_size(1));
        GaborData.step_size = -GaborData.step_size;
    end
end

if isequal(GaborData.stair_fn, @Staircase.noise)
    if GaborData.step_size(1) < 0
        warning('Changing sign of noise step_size from %d to %d', GaborData.step_size(1), -GaborData.step_size(1));
        GaborData.step_size = -GaborData.step_size;
    end
    if ~(iseffectiveinteger(GaborData.step_size(1)) && iseffectiveinteger(GaborData.min_step_size))
        error('In Noise condition, step_size and min_step_size act on indices and must be integers');
    end
end

if isequal(GaborData.stair_fn, @Staircase.contrast)
    if GaborData.step_size(1) < 0
        error('Contrast staircase is multiplicative; step size of %f doesn''t make sense', GaborData.step_size(1));
    elseif GaborData.step_size(1) < 1
        warning('Chaning contrast step_size < 1 to 1/step_size');
        GaborData.step_size(1) = 1 / GaborData.step_size(1);
    end
end

if GaborData.ratio(1) > 1 || GaborData.ratio(1) < 0
    error('Ratio should be between 0.5 and 1');
end

if ~isempty(GaborData.model_observer) && ~any(strcmpi(GaborData.model_observer, {'ideal', 'oracle', 'bernoulli'}))
    warning('%s is not a known model observer', GaborData.model_observer);
end

if strcmp(GaborData.model_observer, 'bernoulli')
    GaborData.sigmoid_slope = get_arg('sigmoid_slope', 20);
    default_pk = ones(1, GaborData.number_of_images) / GaborData.number_of_images;
    GaborData.model_pk = get_arg('model_pk', default_pk);
    if sum(GaborData.model_pk) ~= 1
        warning('Recommended that GaborData.model_pk sum to 1');
    end
end

if ~isempty(varargin)
    warning('Unkown arguments given to newGaborParams');
end

disp(GaborData);


    function TF = iseffectiveinteger(v)
        TF = (v == floor(v));
    end

end
