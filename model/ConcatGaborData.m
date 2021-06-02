function GaborData = ConcatGaborData(Data1, Data2)
%CONCATPRELIMGABOR Concatenate experiment results.

GaborData = Data1;

GaborData.current_trial = Data1.current_trial + Data2.current_trial;
if size(Data1.contrast,2)==Data1.number_of_images
    Data1.contrast = Data1.contrast(:,1)';
end
GaborData.contrast = [Data1.contrast, Data2.contrast];
% GaborData.contrast_per_frame = [Data1.contrast_per_frame; Data2.contrast_per_frame];
GaborData.ratio = [Data1.ratio, Data2.ratio];
GaborData.noise = [Data1.noise, Data2.noise];
GaborData.step_size = [Data1.step_size, Data2.step_size];

GaborData.iid = [Data1.iid, Data2.iid];
GaborData.seed = [Data1.seed, Data2.seed];
if isfield(GaborData, 'checksum'), GaborData.checksum = [Data1.checksum, Data2.checksum]; end
GaborData.streak = [Data1.streak, Data2.streak];
GaborData.reversal_counter = [Data1.reversal_counter, Data2.reversal_counter ];
GaborData.correct_answer = [Data1.correct_answer, Data2.correct_answer];
GaborData.ideal_answer = [Data1.ideal_answer, Data2.ideal_answer];
GaborData.reaction_time = [Data1.reaction_time, Data2.reaction_time];
GaborData.choice = [Data1.choice, Data2.choice];
GaborData.accuracy = [Data1.accuracy, Data2.accuracy];
GaborData.frame_categories = horzcat(Data1.frame_categories, Data2.frame_categories);
GaborData.ideal_frame_signals = horzcat(Data1.ideal_frame_signals, Data2.ideal_frame_signals);
GaborData.image_array_chosen_index = horzcat(Data1.image_array_chosen_index, Data2.image_array_chosen_index);
GaborData.image_array_notchosen_index = horzcat(Data1.image_array_notchosen_index, Data2.image_array_notchosen_index);
GaborData.notchosen_ideal_frame_signals = horzcat(Data1.notchosen_ideal_frame_signals, Data2.notchosen_ideal_frame_signals);
GaborData.chosen_ideal_frame_signals = horzcat(Data1.chosen_ideal_frame_signals, Data2.chosen_ideal_frame_signals);
GaborData.stim_locations = horzcat(Data1.stim_locations, Data2.stim_locations);
GaborData.eye_tracker_points = horzcat(Data1.eye_tracker_points, Data2.eye_tracker_points);
GaborData.gaze_choice = horzcat(Data1.gaze_choice, Data2.gaze_choice);      
GaborData.total_number_of_images = horzcat(Data1.total_number_of_images,Data2.total_number_of_images);
% 
GaborData.true_ratio = [Data1.true_ratio, Data2.true_ratio];
GaborData.sign_noise = [Data1.sign_noise, Data2.sign_noise];
GaborData.sign_contrast = [Data1.sign_contrast, Data2.sign_contrast];
GaborData.num_peri = [Data1.num_peri(1:Data1.current_trial), Data2.num_peri(1:Data2.current_trial)];
% 
GaborData.eye_tracker_points = horzcat(Data1.eye_tracker_points, Data2.eye_tracker_points);
end
