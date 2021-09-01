function [correct_ideal_signals_chosen_regen,correct_ideal_signals_not_chosen_regen,correct_all_ideal_signals_regen] = regenerate_signals(data,trial)

rng(data.seed(trial), 'twister');
first_frame = rand<= data.ratio(trial);
% disp('regenerated first frame:' )
% disp(first_frame);
match_frames = rand(data.num_peri(trial), data.number_of_images) <= data.ratio(trial);
% disp('regenerated other frames')
% disp(match_frames);
% disp('correct answer in data')
% disp(data.correct_answer(trial))
% disp('categories in data')
% disp([data.frame_categories{trial}])
match_frames = logical([first_frame*true(data.num_peri(trial),1) match_frames]);
if data.correct_answer(trial) == 1
    frame_categories(match_frames) = data.left_category;
    frame_categories(~match_frames) = data.right_category;
else
    frame_categories(~match_frames) = data.left_category;
    frame_categories(match_frames) = data.right_category;
end
% disp('categories regenrated')
% disp(reshape(frame_categories,2,data.number_of_images+1))

rng(data.seed(trial), 'twister');
image_array = bpg.genImages(data.num_peri(trial) * (data.number_of_images + 1), data.stim_size, ...
    data.stim_sp_freq_cpp, data.stim_std_sp_freq_cpp, ...
    frame_categories(:), data.noise(trial), data.annulus);
image_array = uint8(image_array * data.contrast(trial) + 127);
image_array = min(image_array, 255);
image_array = max(image_array, 0);
% checksum = mean(image_array(:));
% checksum_regen = [checksum_regen; checksum];
% checksum_stored = [checksum_stored; data.checksum(trial)];
% disp('computed checksum now:')
% disp([checksum])
% disp('saved checksum: ')
% disp([data.checksum(trial)])
image_array_for_sig = [];
image_array_new = [];
image_array_chosen = [];
image_array_not_chosen = [];


image_array_new = reshape(image_array,data.num_peri(trial),data.number_of_images+1,data.stim_size,data.stim_size);
for pr=2:data.num_peri(trial)
    image_array_new(pr,1,:,:) = image_array_new(1,1,:,:);
end
%to compute ideal answer we store images on screen
image_array_for_sig(1,:,:) = squeeze(image_array_new(1,1,:,:));
image_array_chosen(1,:,:) = image_array_for_sig(1,:,:);

tmp_indx = 2;
ch_indx = 2;
notch_indx = 1;
for fr=2:data.number_of_images+1
    for pr=1:data.num_peri(trial)
        image_array_for_sig(tmp_indx,:,:) = squeeze(image_array_new(pr,fr,:,:));
        tmp_indx = tmp_indx + 1;
        if data.image_array_chosen_index{trial}(fr,pr)==true
            image_array_chosen(ch_indx,:,:) = squeeze(image_array_new(pr,fr,:,:));
            ch_indx = ch_indx + 1;
        else
            image_array_not_chosen(notch_indx,:,:) = squeeze(image_array_new(pr,fr,:,:));
            notch_indx = notch_indx + 1;
            
        end
    end
end
correct_all_ideal_signals_regen = ...
    bpg.getSignal(double(image_array_for_sig) - 127, data.left_category, max(data.noise(trial), .04)) - ...
    bpg.getSignal(double(image_array_for_sig) - 127, data.right_category, max(data.noise(trial), .04));
correct_ideal_signals_chosen_regen =  ...
    bpg.getSignal(double(image_array_chosen) - 127, data.left_category, max(data.noise(trial), .04)) - ...
    bpg.getSignal(double(image_array_chosen) - 127, data.right_category, max(data.noise(trial), .04));
correct_ideal_signals_not_chosen_regen = ...
    bpg.getSignal(double(image_array_not_chosen) - 127, data.left_category, max(data.noise(trial), .04)) - ...
    bpg.getSignal(double(image_array_not_chosen) - 127, data.right_category, max(data.noise(trial), .04));
end