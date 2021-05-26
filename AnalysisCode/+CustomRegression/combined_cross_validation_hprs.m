function [best_hprs] = combined_cross_validation_hprs(subjectID, expt_type, hpr_ridge, hpr_ar1, hpr_curvature, standardize, folds)

% initialize variables
[ParentFolderPath] = fileparts(pwd);
datadir = fullfile(ParentFolderPath, '/RawData');
[num_sub,~] = size(subjectID);

for sub =1:num_sub
    signal_all_actual = [];
    signal_all = [];
    choice_raw = [];
    
    % load data
    [data,~] = LoadAllSubjectData(subjectID{sub},expt_type,datadir);
    frames = data.number_of_images;
    num_frames = 2*frames + 1;
    num_peri = 2;
    % we store the signals w.r.t the number of elements on the periphery
    for k = 1:data.current_trial-1
        if data.num_peri(k)==2
            [chosen_ideal_frame_signals,notchosen_ideal_frame_signals,~] = regenerate_signals(data,k);
            signal_all_actual = [signal_all_actual; chosen_ideal_frame_signals' notchosen_ideal_frame_signals'];
            choice_raw = [choice_raw data.choice(k)];
        end
    end
    signal_all = sign(signal_all_actual);
    disp(['Searching hyperparameters for Subject ' num2str(sub)]);
    [~, log_likelihoods(sub,:,:,:,:)] = CustomRegression.xValidatePK_with_lapseSabya(signal_all, choice_raw, hpr_ridge, hpr_ar1, hpr_curvature, standardize, folds);
    log_likelihood_summed(sub,:,:,:) = mean(log_likelihoods(sub,:,:,:,:),5);
    disp(['Done searching for ' num2str(sub) '/' num2str(num_sub) ' subjects...']);
end
sz = size(log_likelihoods);
% avg_ll =  max(mean(log_likelihood_summed,1));
avg_ll =  mean(log_likelihood_summed,1);
[~, imax] = max(avg_ll(:));
[iRidge, iAR1, iCurve] = ind2sub(sz(2:4), imax);
% Err on the side of less regularization by choosing smoothing that is one order of magnitude less than the best.
% iRidge = max(iRidge-1, 1);
% iAR1 = max(iAR1-1, 1);
% iCurve = max(iCurve-1, 1);
best_hprs = [hpr_ridge(iRidge), hpr_ar1(iAR1), hpr_curvature(iCurve)];

disp('Searching of hyperparameters complete!!');
disp (['Best hyperparameters across subjects is: ' num2str(best_hprs)]);

end