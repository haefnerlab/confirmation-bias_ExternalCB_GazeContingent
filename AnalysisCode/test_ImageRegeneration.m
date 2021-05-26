clear all; close all;clc;
subjects = {'GCV3-subject01'};%'GCV3-subject02';'GCV3-subject03';'GCV3-subject07';'GCV3-subject08';'GCV3-subject11';'GCV3-subject13'};
checksum_regen = [];
checksum_stored = [];
faulty_signal_regen = [];
faulty_signal_stored = [];
correct_signal_regen = [];
faulty_chosen_signal_regen =[];
faulty_chosen_signal_stored = [];
correct_chosen_signal_regen = [];
faulty_not_chosen_signal_regen =[];
faulty_not_chosen_signal_stored = [];
correct_not_chosen_signal_regen = [];
correct_signal_regen_notdoubled = [];
correct_chosen_signal_regen_notdoubled = [];
correct_not_chosen_signal_regen_notdoubled = [];
temp1 = [];
temp2 = [];
[num_sub,~] = size(subjects);

for sub=1:num_sub
    [ParentFolderPath] = fileparts(pwd);
    datadir = fullfile(ParentFolderPath, '/RawData');
    
    signal_chosen_raw = [];
    signal_not_chosen_raw = [];
    signal_all = [];
    choice_raw = [];
    accuracy = 0;
    num_peri = 2;
    ideal_answer = [];
    correct_answer = [];
    % load data
    [data,~] = LoadAllSubjectData(subjects{1},3,datadir);
    disp('Data loaded!');
    
    for trial=1:data.current_trial
        if data.num_peri(trial)==2
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
            checksum = mean(image_array(:));
            checksum_regen = [checksum_regen; checksum];
            checksum_stored = [checksum_stored; data.checksum(trial)];
            % disp('computed checksum now:')
            % disp([checksum])
            % disp('saved checksum: ')
            % disp([data.checksum(trial)])
            image_array_for_sig = [];
            image_array_new = [];
            image_array_chosen = [];
            image_array_not_chosen = [];
            
%             image_array_for_sig = uint8(image_array_for_sig);
%             
%             image_array_chosen = uint8(image_array_chosen);
%             
%             image_array_not_chosen = uint8(image_array_not_chosen);
            
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
            
            
            
            
            t_nd = ...
                bpg.getSignal(image_array - 127, data.left_category, max(data.noise(trial), .04)) - ...
                bpg.getSignal(image_array - 127, data.right_category, max(data.noise(trial), .04));
            
            t_d = ...
                bpg.getSignal(double(image_array) - 127, data.left_category, max(data.noise(trial), .04)) - ...
                bpg.getSignal(double(image_array) - 127, data.right_category, max(data.noise(trial), .04));
            
            
            temp_faulty_regen = ...
                bpg.getSignal(image_array_for_sig - 127, data.left_category, max(data.noise(trial), .04)) - ...
                bpg.getSignal(image_array_for_sig - 127, data.right_category, max(data.noise(trial), 4.04));
            
            temp_correct_regen = ...
                bpg.getSignal(double(image_array_for_sig) - 127, data.left_category, max(data.noise(trial), .04)) - ...
                bpg.getSignal(double(image_array_for_sig) - 127, data.right_category, max(data.noise(trial), .04));
            
            temp_correct_regen_notdoubled = ...
                bpg.getSignal(image_array_for_sig - 127, data.left_category, max(data.noise(trial), .04)) - ...
                bpg.getSignal(image_array_for_sig - 127, data.right_category, max(data.noise(trial), .04));
            
            
            temp_stored = data.ideal_frame_signals{trial};
            
            temp_faulty_signal_chosen_regen =  ...
                bpg.getSignal(image_array_chosen - 127, data.left_category, max(data.noise(trial), .04)) - ...
                bpg.getSignal(image_array_chosen - 127, data.right_category, max(data.noise(trial), 4.04));
            
            temp_correct_signal_chosen_regen =  ...
                bpg.getSignal(double(image_array_chosen) - 127, data.left_category, max(data.noise(trial), .04)) - ...
                bpg.getSignal(double(image_array_chosen) - 127, data.right_category, max(data.noise(trial), .04));
            
            temp_correct_signal_chosen_regen_notdoubled =  ...
                bpg.getSignal(image_array_chosen - 127, data.left_category, max(data.noise(trial), .04)) - ...
                bpg.getSignal(image_array_chosen - 127, data.right_category, max(data.noise(trial), .04));
            
            
            temp_faulty_signal_not_chosen_regen = ...
                bpg.getSignal(image_array_not_chosen - 127, data.left_category, max(data.noise(trial), .04)) - ...
                bpg.getSignal(image_array_not_chosen - 127, data.right_category, max(data.noise(trial), 4.04));
            
            temp_correct_signal_not_chosen_regen = ...
                bpg.getSignal(double(image_array_not_chosen) - 127, data.left_category, max(data.noise(trial), .04)) - ...
                bpg.getSignal(double(image_array_not_chosen) - 127, data.right_category, max(data.noise(trial), .04));
            
            temp_correct_signal_not_chosen_regen_notdoubled = ...
                bpg.getSignal(image_array_not_chosen - 127, data.left_category, max(data.noise(trial), .04)) - ...
                bpg.getSignal(image_array_not_chosen - 127, data.right_category, max(data.noise(trial), .04));
            
            
            temp_signal_chosen_stored = data.chosen_ideal_frame_signals{trial};
            
            temp_signal_not_chosen_stored = data.notchosen_ideal_frame_signals{trial};
            
            faulty_signal_regen =[faulty_signal_regen; temp_faulty_regen(:)];
            faulty_signal_stored = [faulty_signal_stored; temp_stored(:)];
            correct_signal_regen = [correct_signal_regen; temp_correct_regen(:)];
            correct_signal_regen_notdoubled = [correct_signal_regen_notdoubled; temp_correct_regen_notdoubled(:)];
            
            faulty_chosen_signal_regen =[faulty_chosen_signal_regen; temp_faulty_signal_chosen_regen(:)];
            faulty_chosen_signal_stored = [faulty_chosen_signal_stored; temp_signal_chosen_stored(:)];
            correct_chosen_signal_regen = [correct_chosen_signal_regen; temp_correct_signal_chosen_regen(:)];
            correct_chosen_signal_regen_notdoubled = [correct_chosen_signal_regen_notdoubled; temp_correct_signal_chosen_regen_notdoubled(:)];
            
            faulty_not_chosen_signal_regen =[faulty_not_chosen_signal_regen; temp_faulty_signal_not_chosen_regen(:)];
            faulty_not_chosen_signal_stored = [faulty_not_chosen_signal_stored; temp_signal_not_chosen_stored(:)];
            correct_not_chosen_signal_regen = [correct_not_chosen_signal_regen; temp_correct_signal_not_chosen_regen(:)];
            correct_not_chosen_signal_regen_notdoubled = [correct_not_chosen_signal_regen_notdoubled; temp_correct_signal_not_chosen_regen_notdoubled(:)];
 
            temp1 = [temp1; t_d];
            temp2 = [temp2; t_nd];
        end
    end
end
%%
figure();
vals1 = linspace(min(min(checksum_regen),min(checksum_stored)),max(max(checksum_regen),max(checksum_stored)),100);
subplot(1,3,1)
scatter(checksum_stored,checksum_regen,'o');
hold on;
plot(vals1,vals1,'k');
xlabel('Saved mean of image array per trial')
ylabel('Recovered mean of image array per trial')
title('Matching mean of image array')

vals2 = linspace(min(min(faulty_signal_regen),min(faulty_signal_stored)),max(max(faulty_signal_regen),max(faulty_signal_stored)),100);
subplot(1,3,2)
scatter(faulty_signal_stored,faulty_signal_regen,'o');
hold on;
plot(vals2,vals2,'k');
xlabel('Saved faulty ideal signals')
ylabel('Recovered faulty ideal signals')
title('Matching faulty signals')

vals3 = linspace(min(min(correct_signal_regen),min(faulty_signal_stored)),max(max(correct_signal_regen),max(faulty_signal_stored)),100);
subplot(1,3,3)
scatter(faulty_signal_stored,correct_signal_regen,'o');
hold on;
plot(vals3,vals3,'k');
xlabel('Saved faulty signals')
ylabel('Recovered correct signals')
title('Comparing faulty stored and correct recovered signals')
%%

figure();
vals1 = linspace(min(min(faulty_chosen_signal_regen),min(faulty_chosen_signal_stored)),max(max(faulty_chosen_signal_regen),max(faulty_chosen_signal_stored)),100);
subplot(2,2,1)
scatter(faulty_chosen_signal_stored,faulty_chosen_signal_regen,'o');
hold on;
plot(vals1,vals1,'k');
xlabel('Saved faulty chosen signals')
ylabel('Recovered faulty chosen signals')
title('Matching faulty chosen signals')

vals2 = linspace(min(min(faulty_not_chosen_signal_regen),min(faulty_not_chosen_signal_stored)),max(max(faulty_not_chosen_signal_regen),max(faulty_not_chosen_signal_stored)),100);
subplot(2,2,2)
scatter(faulty_not_chosen_signal_stored,faulty_not_chosen_signal_regen,'o');
hold on;
plot(vals2,vals2,'k');
xlabel('Saved faulty not chosen signals')
ylabel('Recovered faulty not chosen signals')
title('Matching faulty not chosen signals')

vals3 = linspace(min(min(correct_chosen_signal_regen),min(faulty_chosen_signal_stored)),max(max(correct_chosen_signal_regen),max(faulty_chosen_signal_stored)),100);
subplot(2,2,3)
scatter(faulty_chosen_signal_stored,correct_chosen_signal_regen,'o');
hold on;
plot(vals3,vals3,'k');
xlabel('Saved faulty chosen signals')
ylabel('Recovered correct chosen signals')
title('Comparing faulty stored and correct recovered chosen signals')

vals4 = linspace(min(min(correct_not_chosen_signal_regen),min(faulty_not_chosen_signal_stored)),max(max(correct_not_chosen_signal_regen),max(faulty_not_chosen_signal_stored)),100);
subplot(2,2,4)
scatter(faulty_not_chosen_signal_stored,correct_not_chosen_signal_regen,'o');
hold on;
plot(vals4,vals4,'k');
xlabel('Saved faulty not chosen signals')
ylabel('Recovered correct not chosen signals')
title('Comparing faulty stored and correct recovered not chosen signals')

%%
figure();
vals1 = linspace(min(min(correct_chosen_signal_regen),min(correct_chosen_signal_regen_notdoubled)),max(max(correct_chosen_signal_regen),max(correct_chosen_signal_regen_notdoubled)),100);
subplot(1,4,1)
scatter(correct_chosen_signal_regen_notdoubled,correct_chosen_signal_regen,'o');
hold on;
plot(vals1,vals1,'k');
xlabel('Not doubled chosen signals')
ylabel('Doubled correct chosen signals')
title('Comparing doubled and not doubled chosen signals')

vals2 = linspace(min(min(correct_not_chosen_signal_regen),min(correct_not_chosen_signal_regen_notdoubled)),max(max(correct_not_chosen_signal_regen),max(correct_not_chosen_signal_regen_notdoubled)),100);
subplot(1,4,2)
scatter(correct_not_chosen_signal_regen_notdoubled,correct_not_chosen_signal_regen,'o');
hold on;
plot(vals2,vals2,'k');
xlabel('Not doubled not chosen signals')
ylabel('Doubled correct not chosen signals')
title('Comparing doubled and not doubled not chosen signals')

vals3 = linspace(min(min(correct_signal_regen),min(correct_signal_regen_notdoubled)),max(max(correct_signal_regen),max(correct_signal_regen_notdoubled)),100);
subplot(1,4,3)
scatter(correct_signal_regen_notdoubled,correct_signal_regen,'o');
hold on;
plot(vals3,vals3,'k');
xlabel('Not doubled all signals')
ylabel('Doubled all signals')
title('Comparing doubled and not doubled all signals')

vals4 = linspace(min(min(temp1),min(temp2)),max(max(temp1),max(temp2)),100);
subplot(1,4,4)
scatter(temp1,temp2,'o');
hold on;
plot(vals4,vals4,'k');
xlabel('Not doubled')
ylabel('Doubled')
title('Comparing doubled and not doubled')
