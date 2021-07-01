function [image_array, frame_categories, checksum] = GaborStimulus_staircase(GaborData, trial)
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


% Choose frames based on whether correct answer this trial is Left or Right
if GaborData.correct_answer(trial) == 1
    if rand<=0.5
        frame_categories(1) = GaborData.left_category;
        frame_categories(2) = GaborData.left_category;
    else
        frame_categories(1) = GaborData.right_category;
        frame_categories(2) = GaborData.right_category;
    end
    
else
    if rand<=0.5
        frame_categories(1) = GaborData.left_category;
        frame_categories(2) = GaborData.right_category;
    else
        frame_categories(1) = GaborData.right_category;
        frame_categories(2) = GaborData.left_category;
    end
end
frame_categories_temp = frame_categories(:);
% Set random seed again to keep match_frames independent of pixel noise.
rng(GaborData.seed(trial), 'twister');
image_array = stim_fcn(2, GaborData.stim_size, ...
    GaborData.stim_sp_freq_cpp, GaborData.stim_std_sp_freq_cpp, ...
    frame_categories_temp, GaborData.noise(trial), GaborData.annulus);

image_array = uint8(image_array * GaborData.contrast(trial) + 127);

image_array = min(image_array, 255);
image_array = max(image_array, 0);

checksum = mean(image_array(:));

if isfield(GaborData, 'checksum') && GaborData.checksum(trial) ~= 0 && GaborData.checksum(trial) ~= checksum
    error('Stimulus reconstruction checksum failed!');
end

end