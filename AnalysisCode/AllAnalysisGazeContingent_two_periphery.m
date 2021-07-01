function [params_boot, sobl, abbl, best_hprs,...
    data, accuracy, not_wt_belief_opp_cases,...
    not_wt_trials_opp] = AllAnalysisGazeContingent_two_periphery(subjectID, expt_type, boot_n, gaps, hpr1,hpr2)

% initialize variables
[ParentFolderPath] = fileparts(pwd);
datadir = fullfile(ParentFolderPath, '/RawData');


signal_chosen_raw = [];
signal_not_chosen_raw = [];
signal_all = [];
choice_raw = [];
accuracy = 0;
num_peri = 2;

% load data
[data,~] = LoadAllSubjectData(subjectID,expt_type,datadir);
disp('Data loaded!');

frames = data.number_of_images;


% we store the signals w.r.t the number of elements on the periphery
for k = 1:data.current_trial-1
    if data.num_peri(k)==2
        signal_chosen_raw = [signal_chosen_raw; data.chosen_ideal_frame_signals{k}];
        signal_not_chosen_raw = [signal_not_chosen_raw; data.notchosen_ideal_frame_signals{k}];
        signal_all = [signal_all; data.chosen_ideal_frame_signals{k}' data.notchosen_ideal_frame_signals{k}'];
        choice_raw = [choice_raw data.choice(k)];
        accuracy = accuracy + data.accuracy(k);
    end
end

% this part computes the hyperparameters for the PK
[best_hprs, ~] = xValidatePK_with_lapse(signal_all, choice_raw, frames, hpr1, 0, hpr2, 0, 10);
disp('best params found')
% this part computes the PK weights for the chosen and the not chosen
% signals, to be later used for computation of accumulated evidence
trials = size(choice_raw, 2);
all_frames = 1 + frames * num_peri;
for j = 1:boot_n
    [signal, choice] = bootstrap(signal_all, choice_raw, trials);
    [sobl(j,:), ~] = LinearPK_with_lapse(signal(:,1:frames), choice, 0);
    disp(j)
    [params_boot(j,:), ~, ~, ~, ~, ~] = PsychophysicalKernelwithlapse(signal, choice, frames, best_hprs(1), 0, best_hprs(3), 0);
    [abbl(j,:), ~, ~] = ExponentialPK_with_lapse(signal(:,1:frames), choice, 0);
end

% signals weighted by PK weights based on final categorical choice
weighted_signal_all = (prctile(params_boot(:,1:all_frames),50)).* signal_all;
weighted_signal_chosen = weighted_signal_all(:,1:frames);
weighted_signal_notchosen = weighted_signal_all(:,frames+1:end-num_peri);
weighted_signal_notchosen = reshape(weighted_signal_notchosen,trials,num_peri-1,frames-1);

% ideal signals not weighted at all
signal_chosen = signal_all(:,1:frames);
signal_notchosen = signal_all(:,frames+1:end-num_peri);
signal_notchosen = reshape(signal_notchosen,trials,num_peri-1,frames-1);

% computing the effects of external CB
[not_wt_belief_opp_cases,not_wt_trials_opp]...
    = externalCBAnalysisGazeContingent_not_weighted_choice(weighted_signal_chosen,weighted_signal_notchosen,signal_chosen,signal_notchosen,params_boot(j,end-1),gaps,num_peri,boot_n);


    function [signals, choices] = bootstrap(signals_raw, choices_raw, trials)
        sample_nums = randsample(trials, trials, true); % random sample with replacement
        signals = [];
        choices = [];
        for i = 1:trials
            trial_num = sample_nums(i);
            signals = [signals; signals_raw(trial_num, :)];
            choices = [choices choices_raw(trial_num)];
        end
    end

end