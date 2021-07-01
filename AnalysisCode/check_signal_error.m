clear all; close all;clc;
subjects = {'GCV3-subject01';'GCV3-subject02';'GCV3-subject03';'GCV3-subject07';'GCV3-subject08';'GCV3-subject11';'GCV3-subject13'};
[num_sub,~] = size(subjects);
figure();
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
    [data,~] = LoadAllSubjectData(subjects{sub},3,datadir);
    disp('Data loaded!');
    
    frames = data.number_of_images;
    
    
%     % we store the signals w.r.t the number of elements on the periphery
%     for k = 1:data.current_trial-1
%         if data.num_peri(k)==2
%             signal_chosen_raw = [signal_chosen_raw; data.chosen_ideal_frame_signals{k}];
%             signal_not_chosen_raw = [signal_not_chosen_raw; data.notchosen_ideal_frame_signals{k}];
%             signal_all = [signal_all; data.chosen_ideal_frame_signals{k}' data.notchosen_ideal_frame_signals{k}'];
%             choice_raw = [choice_raw data.choice(k)];
%             accuracy = accuracy + data.accuracy(k);
%             ideal_answer = [ideal_answer data.ideal_answer(k)];
%             correct_answer =[correct_answer data.correct_answer(k)];
%         end
%     end
%     
%     subplot(3,num_sub,sub)
%     hist(signal_all(:));
%     xlabel('signal')
%     subplot(3,num_sub,num_sub+sub)
%     hist(ideal_answer);
%     xlabel('ideal answer')
%     subplot(3,num_sub,2*num_sub+sub)
%     hist(ideal_answer);
%     xlabel('correct answer');
%     hold on;

      
end
