function data_sub = readSubjectData(subjectID)
expt_type = 3;
[ParentFolderPath] = fileparts(pwd);
datadir = fullfile(ParentFolderPath, '/RawData');

signal_chosen_raw_actual = [];
signal_not_chosen_raw_actual = [];
signal_all_actual = [];
choice_raw = [];

% load data
[data,~] = LoadAllSubjectData(subjectID,expt_type,datadir);
disp('Data for the subject loaded!');
frames = data.number_of_images;

accuracy = [];
accuracy_ideal_foveal_only = [];
accuracy_ideal_all = [];
num_peri = 2;
% we store the signals w.r.t the number of elements on the periphery
for k = 1:data.current_trial-1
    if data.num_peri(k)==2
        [chosen_ideal_frame_signals,notchosen_ideal_frame_signals,~] = regenerate_signals(data,k);
        signal_chosen_raw_actual = [signal_chosen_raw_actual; chosen_ideal_frame_signals'];
        signal_not_chosen_raw_actual = [signal_not_chosen_raw_actual; notchosen_ideal_frame_signals'];
        signal_all_actual = [signal_all_actual; chosen_ideal_frame_signals' notchosen_ideal_frame_signals'];
        choice_raw = [choice_raw data.choice(k)];
        accuracy = [accuracy data.accuracy(k)];
        accuracy_ideal_foveal_only = [accuracy_ideal_foveal_only (1 * (sum(chosen_ideal_frame_signals)>0))==data.ideal_answer(k)];
        accuracy_ideal_all = [accuracy_ideal_all (1 * (sum(chosen_ideal_frame_signals)+sum(notchosen_ideal_frame_signals))>0)==data.ideal_answer(k)];
    end
end

signal_chosen_raw = sign(signal_chosen_raw_actual);
signal_not_chosen_raw = sign(signal_not_chosen_raw_actual);
signal_all = sign(signal_all_actual);

trials = size(choice_raw, 2);
all_frames = 1 + frames * num_peri;

data_sub.frame_signals = signal_all_actual;
data_sub.frame_categories = signal_all;
data_sub.choice = choice_raw;
data_sub.num_trials = trials;
data_sub.accuracy = accuracy; 
end