function [params_boot, sobl, abbl, log_bernoulli,...
    data, accuracy, belief_opp_cases, ...
    trials_opp,belief_opp_cases_raw,...
    mn_sig,mn_sig_norm,mn_sig_ideal_foveal_only,mn_sig_ideal_foveal_only_norm,...
    mn_sig_ideal_all,mn_sig_ideal_all_norm,mn_perf,mn_perf_ideal_foveal_only,mn_perf_ideal_all,...
    tr_seg,err_perf,tr_seg_ideal_foveal_only,err_perf_ideal_foveal_only,tr_seg_ideal_all,err_perf_ideal_all,...
    mn_sig_ideal_all_actual_sig,mn_perf_actual_sig,mn_perf_ideal_foveal_only_actual_sig,mn_perf_ideal_all_actual_sig,...
    tr_seg_actual_sig,err_perf_actual_sig,tr_seg_ideal_foveal_only_actual_sig,err_perf_ideal_foveal_only_actual_sig,...
    tr_seg_ideal_all_actual_sig,err_perf_ideal_all_actual_sig,...
    actual_data, predicted_data, mean_sig_shown, prob_opp_divided,...
    belief_opp_cases_first_half,belief_opp_cases_second_half,accuracy_first_half,accuracy_second_half]...
    = AnalysisGazeContingent_SingleHyprs(subjectID, expt_type, boot_n, gaps, best_hprs, segments, sanity_bins, standardize)

% initialize variables
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

trials_first_half = floor(trials/2);
trials_second_half = trials - trials_first_half;

signal_chosen_raw_first_half = sign(signal_chosen_raw_actual(1:trials_first_half,:));
signal_not_chosen_raw_first_half = sign(signal_not_chosen_raw_actual(1:trials_first_half,:));
signal_all_first_half = sign(signal_all_actual(1:trials_first_half,:));
signal_chosen_raw_second_half = sign(signal_chosen_raw_actual(trials_first_half+1:end,:));
signal_not_chosen_raw_second_half = sign(signal_not_chosen_raw_actual(trials_first_half+1:end,:));
signal_all_second_half = sign(signal_all_actual(trials_first_half+1:end,:));
choice_raw_first_half = choice_raw(1:trials_first_half);
choice_raw_second_half = choice_raw(trials_first_half+1:end);
accuracy_first_half = sum(accuracy(1:trials_first_half))/trials_first_half;
accuracy_second_half = sum(accuracy(trials_first_half+1:end))/trials_second_half;

disp('Computing regression weights...');
for j = 1:boot_n
    [signal, choice] = bootstrap(signal_all, choice_raw, trials);
    [signal_first_half, choice_first_half] = bootstrap(signal_all_first_half, choice_raw_first_half, trials_first_half);
    [signal_second_half, choice_second_half] = bootstrap(signal_all_second_half, choice_raw_second_half, trials_second_half);
    %     [sobl(j,:), ~] = CustomRegression.LinearPK_with_lapse(signal(:,1:frames), choice, 0);
    sobl(j,:) = zeros(1,4);
    if mod(j,100)==0 || j==1
        disp(['Bootstrap number ' num2str(j) '/' num2str(boot_n) ' ...']);
    end
    [params_boot(j,:), ~, ~, ~, ~, ~] = CustomRegression.PsychophysicalKernelwithlapseSabya(signal, choice, best_hprs(1), best_hprs(2), best_hprs(3), standardize);
    [params_boot_first_half(j,:), ~, ~, ~, ~, ~] = CustomRegression.PsychophysicalKernelwithlapseSabya(signal_first_half, choice_first_half, best_hprs(1), best_hprs(2), best_hprs(3), standardize);
    [params_boot_second_half(j,:), ~, ~, ~, ~, ~] = CustomRegression.PsychophysicalKernelwithlapseSabya(signal_second_half, choice_second_half, best_hprs(1), best_hprs(2), best_hprs(3), standardize);
    %     [abbl(j,:), ~, ~] = CustomRegression.ExponentialPK_with_lapse(signal(:,1:frames), choice, 0);
    abbl(j,:) = zeros(1,4);
end

% signals weighted by PK weights based on final categorical choice
temporal_kernel = prctile(params_boot(:, 1:all_frames), 50);
weighted_signal_all = temporal_kernel .* signal_all;
weighted_signal_chosen = weighted_signal_all(:,1:frames);
weighted_signal_notchosen = weighted_signal_all(:,frames+1:end-num_peri);
weighted_signal_notchosen = reshape(weighted_signal_notchosen,trials,num_peri-1,frames-1);
bias =  prctile(params_boot(:, end-1), 50);
lapse = prctile(1e-4+(1-1e-4) * sigmoid(params_boot(:,end)),50);
% lapse = prctile((params_boot(:,end)).^2,50);

temporal_kernel_first_half = prctile(params_boot_first_half(:, 1:all_frames), 50);
weighted_signal_all_first_half = temporal_kernel_first_half .* signal_all_first_half;
weighted_signal_chosen_first_half = weighted_signal_all_first_half(:,1:frames);
weighted_signal_notchosen_first_half = weighted_signal_all_first_half(:,frames+1:end-num_peri);
weighted_signal_notchosen_first_half = reshape(weighted_signal_notchosen_first_half,trials_first_half,num_peri-1,frames-1);
bias_first_half =  prctile(params_boot_first_half(:, end-1), 50);
lapse_first_half = prctile(1e-4+(1-1e-4) * sigmoid(params_boot_first_half(:,end)),50);

temporal_kernel_second_half = prctile(params_boot_second_half(:, 1:all_frames), 50);
weighted_signal_all_second_half = temporal_kernel_second_half .* signal_all_second_half;
weighted_signal_chosen_second_half = weighted_signal_all_second_half(:,1:frames);
weighted_signal_notchosen_second_half = weighted_signal_all_second_half(:,frames+1:end-num_peri);
weighted_signal_notchosen_second_half = reshape(weighted_signal_notchosen_second_half,trials_second_half,num_peri-1,frames-1);
bias_second_half =  prctile(params_boot_second_half(:, end-1), 50);
lapse_second_half = prctile(1e-4+(1-1e-4) * sigmoid(params_boot_second_half(:,end)),50);

disp('Getting log odds...');
[log_bernoulli] = compute_log_odds(signal_all, temporal_kernel, bias);

% ideal signals not weighted at all
signal_chosen = sign(signal_all(:,1:frames));
signal_notchosen = sign(signal_all(:,frames+1:end-num_peri));
signal_notchosen = reshape(signal_notchosen,trials,num_peri-1,frames-1);

signal_chosen_first_half = sign(signal_all_first_half(:,1:frames));
signal_notchosen_first_half = sign(signal_all_first_half(:,frames+1:end-num_peri));
signal_notchosen_first_half = reshape(signal_notchosen_first_half,trials_first_half,num_peri-1,frames-1);

signal_chosen_second_half = sign(signal_all_second_half(:,1:frames));
signal_notchosen_second_half = sign(signal_all_second_half(:,frames+1:end-num_peri));
signal_notchosen_second_half = reshape(signal_notchosen_second_half,trials_second_half,num_peri-1,frames-1);



disp('Computing bias in saccade...');
% computing the effects of external CB
for gp=1:length(gaps)
    %     [not_wt_belief_opp_cases{gp},not_wt_trials_opp{gp}]...
    %         = externalCBAnalysisGazeContingent_not_weighted_choice_bootstrap(weighted_signal_chosen,weighted_signal_notchosen,signal_chosen,signal_notchosen,bias,gaps(gp),num_peri,boot_n);
    [belief_opp_cases{gp},trials_opp{gp}]...
        = externalCBAnalysisGazeContingent(weighted_signal_chosen,weighted_signal_notchosen,signal_chosen,signal_notchosen,bias,gaps(gp),num_peri);
    [belief_opp_cases_raw{gp}] = externalCBAnalysisGazeContingent_raw(signal_chosen,signal_notchosen,gaps(gp),num_peri,boot_n);
end
[prob_opp_divided] = externalCBAnalysisGazeContingent_divided(weighted_signal_chosen,weighted_signal_notchosen,signal_chosen,signal_notchosen,bias,num_peri);
[belief_opp_cases_first_half,~]...
        = externalCBAnalysisGazeContingent(weighted_signal_chosen_first_half,weighted_signal_notchosen_first_half,signal_chosen_first_half,signal_notchosen_first_half,bias_first_half,2,num_peri);
[belief_opp_cases_second_half,~]...
        = externalCBAnalysisGazeContingent(weighted_signal_chosen_second_half,weighted_signal_notchosen_second_half,signal_chosen_second_half,signal_notchosen_second_half,bias_second_half,2,num_peri);
    
disp('Computing performance w.r.t evidence accumulated...');
[mn_sig,mn_sig_norm,mn_sig_ideal_foveal_only,mn_sig_ideal_foveal_only_norm,...
    mn_sig_ideal_all,mn_sig_ideal_all_norm,mn_perf,mn_perf_ideal_foveal_only,mn_perf_ideal_all,...
    tr_seg,err_perf,tr_seg_ideal_foveal_only,err_perf_ideal_foveal_only,tr_seg_ideal_all,err_perf_ideal_all]...
    = compute_performance_two_periphery(temporal_kernel,signal_all,signal_chosen,...
    accuracy,accuracy_ideal_all,accuracy_ideal_foveal_only,all_frames,segments);

disp('Computing performance w.r.t raw signal shown...');
[mn_sig_ideal_all_actual_sig,mn_perf_actual_sig,mn_perf_ideal_foveal_only_actual_sig,mn_perf_ideal_all_actual_sig,...
    tr_seg_actual_sig,err_perf_actual_sig,tr_seg_ideal_foveal_only_actual_sig,err_perf_ideal_foveal_only_actual_sig,...
    tr_seg_ideal_all_actual_sig,err_perf_ideal_all_actual_sig]...
    = compute_performance_two_periphery_wrt_actual_signal(signal_all,...
    accuracy,accuracy_ideal_all,accuracy_ideal_foveal_only,segments);

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
    function [logits] = compute_log_odds(data,weights,bias_computed)
        logits = data * weights(:) + bias_computed;
    end

end