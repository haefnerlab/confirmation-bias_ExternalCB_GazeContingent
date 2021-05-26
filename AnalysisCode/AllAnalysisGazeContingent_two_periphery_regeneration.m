function [params_boot, sobl, abbl, best_hprs,...
    data, accuracy, not_wt_belief_opp_cases, ...
    not_wt_trials_opp,belief_opp_cases_raw,...
    mn_sig,mn_sig_norm,mn_sig_ideal_foveal_only,mn_sig_ideal_foveal_only_norm,...
    mn_sig_ideal_all,mn_sig_ideal_all_norm,mn_perf,mn_perf_ideal_foveal_only,mn_perf_ideal_all,...
    tr_seg,err_perf,tr_seg_ideal_foveal_only,err_perf_ideal_foveal_only,tr_seg_ideal_all,err_perf_ideal_all,...
    actual_data, predicted_data, mean_sig_shown]...
    = AllAnalysisGazeContingent_two_periphery_regeneration(subjectID, expt_type, boot_n, gaps, hpr1, hpr2, segments, sanity_bins)

% initialize variables
[ParentFolderPath] = fileparts(pwd);
datadir = fullfile(ParentFolderPath, '/RawData');

signal_chosen_raw = [];
signal_not_chosen_raw = [];
signal_all = [];
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
        signal_chosen_raw = [signal_chosen_raw; chosen_ideal_frame_signals];
        signal_not_chosen_raw = [signal_not_chosen_raw; notchosen_ideal_frame_signals];
        signal_all = [signal_all; chosen_ideal_frame_signals' notchosen_ideal_frame_signals'];
        choice_raw = [choice_raw data.choice(k)];
        accuracy = [accuracy data.accuracy(k)];
        accuracy_ideal_foveal_only = [accuracy_ideal_foveal_only (1 * (sum(chosen_ideal_frame_signals)>0))==data.ideal_answer(k)];
        accuracy_ideal_all = [accuracy_ideal_all (1 * (sum(chosen_ideal_frame_signals)+sum(notchosen_ideal_frame_signals))>0)==data.ideal_answer(k)];
        
    end
end
disp('Searching hyperparameters...');
% this part computes the hyperparameters for the PK
[best_hprs, ~] = CustomRegression.xValidatePK_with_lapseSabya(signal_all, choice_raw, frames, hpr1, 0, hpr2, 0, 10);
% [best_hprs, ~] = CustomRegression.xValidatePK_with_lapse(signal_all, choice_raw, frames, hpr1, 0, hpr2, 0, 10, opt_type);
disp(['Best hyperparameters found: ' num2str(best_hprs)]);
% this part computes the PK weights for the chosen and the not chosen
% signals, to be later used for computation of accumulated evidence
trials = size(choice_raw, 2);
all_frames = 1 + frames * num_peri;
disp('Computing regression weights...');
for j = 1:boot_n
    [signal, choice] = bootstrap(signal_all, choice_raw, trials);
%     [sobl(j,:), ~] = CustomRegression.LinearPK_with_lapse(signal(:,1:frames), choice, 0);
    sobl(j,:) = zeros(1,4);
    if mod(j,200)==0 || j==1
        disp(['Bootstrap number ' num2str(j) '/' num2str(boot_n) ' ...']);
    end
    [params_boot(j,:), ~, ~, ~, ~, ~] = CustomRegression.PsychophysicalKernelwithlapseSabya(signal, choice, frames, best_hprs(1), 0, best_hprs(3), 0);
%     [params_boot(j,:), ~, ~, ~, ~, ~] = CustomRegression.PsychophysicalKernelwithlapse(signal, choice, frames, best_hprs(1), 0, best_hprs(3), 0, opt_type);
    %     [abbl(j,:), ~, ~] = CustomRegression.ExponentialPK_with_lapse(signal(:,1:frames), choice, 0);
    abbl(j,:) = zeros(1,4);
end

% signals weighted by PK weights based on final categorical choice
weighted_signal_all = prctile(params_boot(:, 1:all_frames), 50).* signal_all;
weighted_signal_chosen = weighted_signal_all(:,1:frames);
weighted_signal_notchosen = weighted_signal_all(:,frames+1:end-num_peri);
bias =  prctile(params_boot(:, end-1), 50);
% lapse = prctile(1e-4+(1-1e-4)*sigmoid(params_boot(:,end)),50);
lapse = prctile((params_boot(:,end)).^2,50);

% shuffled weighted signal for baseline comparison
weighted_signal_notchosen = reshape(weighted_signal_notchosen,trials,num_peri-1,frames-1);

% ideal signals not weighted at all
signal_chosen = signal_all(:,1:frames);
signal_notchosen = signal_all(:,frames+1:end-num_peri);
signal_notchosen = reshape(signal_notchosen,trials,num_peri-1,frames-1);

disp('Computing bias in saccade...');
% computing the effects of external CB
for gp=1:length(gaps)
    [not_wt_belief_opp_cases{gp},not_wt_trials_opp{gp}]...
        = externalCBAnalysisGazeContingent_not_weighted_choice(weighted_signal_chosen,weighted_signal_notchosen,signal_chosen,signal_notchosen,params_boot(j,end-1),gaps(gp),num_peri,boot_n);
    [belief_opp_cases_raw{gp}] = externalCBAnalysisGazeContingent_raw(signal_chosen,signal_notchosen,gaps(gp), num_peri, boot_n);
end

disp('Computing performance w.r.t evidence accumulated...');
temporal_kernel = prctile(params_boot(:, 1:all_frames), 50);
[mn_sig,mn_sig_norm,mn_sig_ideal_foveal_only,mn_sig_ideal_foveal_only_norm,...
    mn_sig_ideal_all,mn_sig_ideal_all_norm,mn_perf,mn_perf_ideal_foveal_only,mn_perf_ideal_all,...
    tr_seg,err_perf,tr_seg_ideal_foveal_only,err_perf_ideal_foveal_only,tr_seg_ideal_all,err_perf_ideal_all]...
    = compute_performance_two_periphery(temporal_kernel,signal_all,signal_chosen,...
    accuracy,accuracy_ideal_all,accuracy_ideal_foveal_only,all_frames,segments);

disp('Computing prediction vs true to compare...');
[actual_data, predicted_data, mean_sig_shown] = check_weighted_prediction_vs_true_prediction(choice_raw, signal_all, weighted_signal_all, bias, lapse, sanity_bins);
    
    
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