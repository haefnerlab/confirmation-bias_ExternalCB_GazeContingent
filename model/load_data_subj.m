function [frame_signals_ord,frame_signals,final_choice,category_trial]=load_data_subj(subid)

datadir = fullfile(pwd, '/RawData');
subjectID=['GCV3-subject0',num2str(subid)];
expt_type=3;
signal_chosen_raw = [];
signal_not_chosen_raw = [];
signal_all = [];
choice_raw = [];
choice_corr = [];
accuracy1 = 0;
frame_signals=[];
num_peri = 2;

% load data
[data,~] = LoadAllSubjectData(subjectID,expt_type,datadir);
disp('Data loaded!');

frames = data.number_of_images;


% we store the signals w.r.t the number of elements on the periphery
for k = 1:data.current_trial-1
    if data.num_peri(k)==2
        [chosen_ideal_frame_signals,notchosen_ideal_frame_signals,allsigs] = regenerate_signals(data,k);
        signal_chosen_raw = [signal_chosen_raw; chosen_ideal_frame_signals];
        frame_signals = [frame_signals; allsigs(:)'];
        signal_not_chosen_raw = [signal_not_chosen_raw; notchosen_ideal_frame_signals];
        signal_all = [signal_all; chosen_ideal_frame_signals' notchosen_ideal_frame_signals'];
        choice_raw = [choice_raw data.choice(k)];
        choice_corr = [choice_corr data.correct_answer(k)];
        accuracy1 = accuracy1 + data.accuracy(k);
    end
end

frame_signals_ord=signal_all;
final_choice=sign(choice_raw(:)-0.5);
category_trial=sign(choice_corr(:)-0.5);
end