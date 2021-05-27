function [image_array, frame_categories, checksum] = GaborStimulus(GaborData, trial)
%GABORSTIMULUS(GaborData, trial) create (or recreate) stimulus frames based
%on parameters in GaborData and the seed, contrast, ratio, and noise on the
%given trial. If 'GaborData.iid(trial)' is true, each frame's category is
%drawn iid based on the 'ratio' parameter. Otherwise, exactly
%round(ratio*num_images) frames will match the 'true' category.
%
%This function makes no modifications to GaborData.

% Set RNG state to recreate stimulus for this trail.
rng(GaborData.seed(trial), 'twister');

if isfield(GaborData, 'flag_use_old_stimulus_code') && GaborData.flag_use_old_stimulus_code
    stim_fcn = @bpg.genImagesOld;
else
    stim_fcn = @bpg.genImages;
end

if ~isfield(GaborData, 'iid') || GaborData.iid(trial)
    % Randomly set each frame to match (or mismatch) the correct choice
    % for this trail, using the current 'ratio' to decide.
    first_frame = rand<= GaborData.ratio(trial);
    match_frames = rand(GaborData.num_peri(trial), GaborData.number_of_images) <= GaborData.ratio(trial);
else
    % Randomly permute whether each frame matches the true category, with
    % 'ratio' percent of them matching.
    first_frame = rand<= GaborData.ratio(trial);
    n_match = round(GaborData.ratio(trial) * (GaborData.num_peri(trial) * GaborData.number_of_images));
    match_frames = [true(1, n_match) false(1, (GaborData.num_peri(trial) * GaborData.number_of_images) - n_match)];
    match_frames = Shuffle(match_frames);
    match_frames = reshape(match_frames,GaborData.num_peri(trial),GaborData.number_of_images);
    
end
match_frames = logical([first_frame*true(GaborData.num_peri(trial),1) match_frames]);

frame_categories = zeros(size(match_frames));

% Choose frames based on whether correct answer this trial is Left or Right
if GaborData.correct_answer(trial) == 1
    frame_categories(match_frames) = GaborData.left_category;
    frame_categories(~match_frames) = GaborData.right_category;
else
    frame_categories(~match_frames) = GaborData.left_category;
    frame_categories(match_frames) = GaborData.right_category;
end

% Set random seed again to keep match_frames independent of pixel noise.
rng(GaborData.seed(trial), 'twister');
image_array = stim_fcn(GaborData.num_peri(trial) * (GaborData.number_of_images + 1), GaborData.stim_size, ...
    GaborData.stim_sp_freq_cpp, GaborData.stim_std_sp_freq_cpp, ...
    frame_categories(:), GaborData.noise(trial), GaborData.annulus);

image_array = uint8(image_array * GaborData.contrast(trial) + 127);

image_array = min(image_array, 255);
image_array = max(image_array, 0);

checksum = mean(image_array(:));

if isfield(GaborData, 'checksum') && GaborData.checksum(trial) ~= 0 && GaborData.checksum(trial) ~= checksum
    error('Stimulus reconstruction checksum failed!');
end

end