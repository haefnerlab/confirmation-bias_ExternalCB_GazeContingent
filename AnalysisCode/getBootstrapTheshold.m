function thresh = getBootstrapTheshold(data,performance,boots)
signal_raw = [];
choice_raw = [];
sign_noise_raw = [];
for k = 1:size(data.ideal_frame_signals, 1)
    signal_raw = [signal_raw; data.ideal_frame_signals(k, :)];
    choice_raw = [choice_raw data.choice(k)];
    sign_noise_raw = [sign_noise_raw data.sign_noise(k)];   
end
trials = size(choice_raw, 2);
for i=1:boots
    [signal, choice, sign_noise] = bootstrap(signal_raw, choice_raw,sign_noise_raw ,trials);
    data_temp.choice = choice;
    data_temp.sign_noise = sign_noise;
    data_temp.ideal_frame_signals = signal;
    [pm_fit,~,~,~] = GaborPsychometric(data_temp, -2);
    thresh(i) = getThreshold(pm_fit,performance);
    
end
    function [signal, choice, sign_noise] = bootstrap(signal_raw, choice_raw,sign_noise_raw,trials)
        sample_nums = randsample(trials, trials, true); % random sample with replacement
        signal = [];
        choice = [];
        sign_noise = [];
        for j = 1:trials
            trial_num = sample_nums(j);
            signal = [signal; signal_raw(trial_num, :)];
            choice = [choice choice_raw(trial_num)];
            sign_noise = [sign_noise sign_noise_raw(trial_num)];
        end
    end
end