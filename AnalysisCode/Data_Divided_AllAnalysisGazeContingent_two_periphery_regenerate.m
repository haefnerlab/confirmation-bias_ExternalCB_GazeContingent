function [belief_opp_cases, trials_opp, params_per_session,num_sessions, num_trials_per_session] = Data_Divided_AllAnalysisGazeContingent_two_periphery_regenerate(subjectID, phase, boot_n, gaps, hpr1, hpr2)

% initialize variables
[ParentFolderPath] = fileparts(pwd);
datadir = fullfile(ParentFolderPath, '/RawData');

if phase == 0
    expt_type = 'Contrast';
elseif phase == 1
    expt_type = 'Ratio';
elseif phase == 2
    expt_type = 'Noise';
elseif phase == 3
    expt_type = 'NoStaircase';
else
    error('Expected phase 0 for Contrast or 1 for Ratio or 2 for Noise');
end
num_sessions = 0;
num_trials_per_session = [];
files = dir(fullfile(datadir, '*.mat'));
for i=1:length(files)
    if startsWith(files(i).name, subjectID)
        if endsWith(files(i).name, ['GaborData' expt_type 'Quit.mat'])
            contents = load(fullfile(datadir, files(i).name));
            if contents.GaborData.current_trial < 10
                continue;
            end
            data = TruncateQuitDataGabor(contents.GaborData);
            num_sessions = num_sessions + 1;
            signal_chosen_raw = [];
            signal_not_chosen_raw = [];
            signal_all = [];
            choice_raw = [];
            accuracy = 0;
            num_peri = 2;
            frames = data.number_of_images;
            % we store the signals w.r.t the number of elements on the periphery
            for k = 1:data.current_trial-1
                if data.num_peri(k)==2
                    [chosen_ideal_frame_signals,notchosen_ideal_frame_signals,~] = regenerate_signals(data,k);
                    signal_chosen_raw = [signal_chosen_raw; chosen_ideal_frame_signals];
                    signal_not_chosen_raw = [signal_not_chosen_raw; notchosen_ideal_frame_signals];
                    signal_all = [signal_all; chosen_ideal_frame_signals' notchosen_ideal_frame_signals'];
                    choice_raw = [choice_raw data.choice(k)];
                    accuracy = accuracy + data.accuracy(k);
                end
            end
            
            % this part computes the PK weights for the chosen and the not chosen
            % signals, to be later used for computation of accumulated evidence
            [best_hprs, ~] = xValidatePK_with_lapse(signal_all, choice_raw, frames, hpr1, 0, hpr2, 0, 10);
            disp('best params found')
            % this part computes the PK weights for the chosen and the not chosen
            % signals, to be later used for computation of accumulated evidence
            trials = size(choice_raw, 2);
            num_trials_per_session(end+1) = trials;
            all_frames = 1 + frames * num_peri;
            for j = 1:boot_n
                [signal, choice] = bootstrap(signal_all, choice_raw, trials);
                disp(j)
                [params_boot(j,:), ~, ~, ~, ~, ~] = PsychophysicalKernelwithlapse(signal, choice, frames, best_hprs(1), 0, best_hprs(3), 0);
            end
            params_per_session(num_sessions,:,:) = params_boot;
            % signals weighted by PK weights based on final categorical choice
            for a=1:20
                [sig_wts(a,:), ~, ~, ~, ~, ~] = PsychophysicalKernelwithlapse(signal_all, choice_raw, frames, best_hprs(1), 0, best_hprs(3), 0);
            end
            % signals weighted by PK weights based on final categorical choice
            weighted_signal_all = prctile(sig_wts(:, 1:all_frames), 50).* signal_all;
            weighted_signal_chosen = weighted_signal_all(:,1:frames);
            weighted_signal_notchosen = weighted_signal_all(:,frames+1:end-num_peri);
            weighted_signal_notchosen = reshape(weighted_signal_notchosen,trials,num_peri-1,frames-1);
            
            % ideal signals not weighted at all
            signal_chosen = signal_all(:,1:frames);
            signal_notchosen = signal_all(:,frames+1:end-num_peri);
            signal_notchosen = reshape(signal_notchosen,trials,num_peri-1,frames-1);
            
            % computing the effects of external CB
            [belief_opp_cases(num_sessions,:,:),trials_opp(num_sessions)]...
                = externalCBAnalysisGazeContingent_not_weighted_choice(weighted_signal_chosen,weighted_signal_notchosen,signal_chosen,signal_notchosen,params_boot(j,end-1),gaps,num_peri,boot_n);
            
        else
            continue;
        end
    end
end
    function [signals, choices] = bootstrap(signals_raw, choices_raw, trials)
        sample_nums = randsample(trials, trials, true); % random sample with replacement
        signals = [];
        choices = [];
        for tr = 1:trials
            trial_num = sample_nums(tr);
            signals = [signals; signals_raw(trial_num, :)];
            choices = [choices choices_raw(trial_num)];
        end
    end
end