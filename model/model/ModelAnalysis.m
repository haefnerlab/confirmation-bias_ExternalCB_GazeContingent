clear all; close all;clc;
simulated_data = load('summary_simulation1.mat');
simulated_data = simulated_data.summary;
num_sub = 1;%size(simulated_data,2);
num_peri = 2;
frames = 4;
num_images = num_peri * frames +1;

hpr1 = 0.0;
hpr2 = logspace(-1, 5, 7);
points = 3;
boot_n = 10;

not_wt_belief_extCB_opp = zeros(num_sub,points);
not_wt_belief_err_opp = zeros(num_sub,points);
not_wt_lobelief_err_opp = zeros(num_sub,points);
not_wt_hibelief_err_opp = zeros(num_sub,points);
not_wt_probability_chose_in_favor_opp = zeros(num_sub,points);

figure();
for sub=1:num_sub
    disp(sub)
    GaborData = newGaborData();
    disp(sub);
    sub_data = simulated_data{1,sub};
    num_trials = size(sub_data.cat_trials,1);
    accuracy = 0;
    for tr=1:num_trials
        disp(tr);
        if sub_data.final_choice(tr)==sub_data.cat_trials(tr)
            accuracy = accuracy + 1;
        end
        frame_categories = [];
        frame_categories(1) = sub_data.cat_first_frame(tr);
        for t=1:3
            frame_categories(t+1) = sub_data.cat_frames(tr,sub_data.chosen_loc(tr),t);
            if sub_data.chosen_loc(tr,t)==1
                frame_categories(frames+t) = sub_data.cat_frames(tr,2,t);
            else
                frame_categories(frames+t) = sub_data.cat_frames(tr,1,t);
            end
        end
        frame_categories(end+1:end+2) = sub_data.cat_frames(tr,:,end);
        
        im = bpg.genImages(length(frame_categories), GaborData.stim_size, GaborData.stim_sp_freq_cpp,...
            GaborData.stim_std_sp_freq_cpp, frame_categories, GaborData.noise(1), GaborData.annulus);
        
        im = uint8(im * GaborData.contrast(1) + 127);
        im = min(im, 255);
        im = max(im, 0);
        
        signal_all(tr,:) = bpg.getSignal(double(im) - 127, GaborData.left_category, max(GaborData.noise(1), .04)) - ...
            bpg.getSignal(double(im) - 127, GaborData.right_category, max(GaborData.noise(1), .04));
%         signal_all(tr,:) = bpg.getSignal(double(im) - 127, GaborData.left_category, max(GaborData.noise(1), .04), GaborData.stim_sp_freq_cpp, GaborData.stim_std_sp_freq_cpp) - ...
%             bpg.getSignal(double(im) - 127, GaborData.right_category, max(GaborData.noise(1), .04), GaborData.stim_sp_freq_cpp, GaborData.stim_std_sp_freq_cpp);
%         signal_all_faulty(tr,:) = bpg.getSignal(uint8(im) - 127, GaborData.left_category, max(GaborData.noise(1), .04)) - ...
%             bpg.getSignal(uint8(im) - 127, GaborData.right_category, max(GaborData.noise(1), .04));
        
        if sub_data.final_choice(tr)==0
            choice_raw(tr) = 0;
        else
            choice_raw(tr) = 1;
        end
    end
    disp('Data obtained');
    subplot(1,num_sub,sub)
    hist(signal_all(:));
    hold on;
    %     % this part computes the hyperparameters for the PK
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
        %         [abbl(j,:), ~, ~] = ExponentialPK_with_lapse(signal(:,1:frames), choice, 0);
    end
    [sig_wts, ~, ~, ~, ~, ~] = PsychophysicalKernelwithlapse(signal_all, choice_raw, frames, best_hprs(1), 0, best_hprs(3), 0);
    
    % signals weighted by PK weights based on final categorical choice
    weighted_signal_all = (sig_wts(1:all_frames))'.* signal_all;
    
    weighted_signal_chosen = weighted_signal_all(:,1:frames);
    weighted_signal_notchosen = weighted_signal_all(:,frames+1:end-num_peri);
    weighted_signal_notchosen = reshape(weighted_signal_notchosen,trials,num_peri-1,frames-1);
    
    
    % ideal signals not weighted at all
    signal_chosen = signal_all(:,1:frames);
    signal_notchosen = signal_all(:,frames+1:end-num_peri);
    signal_notchosen = reshape(signal_notchosen,trials,num_peri-1,frames-1);
    
    % computing the effects of external CB
    [not_wt_belief_opp_cases,not_wt_trials_opp]...
        = externalCBAnalysisGazeContingent_not_weighted_choice(weighted_signal_chosen,weighted_signal_notchosen,signal_chosen,signal_notchosen,params_boot(j,end-1),points+1,num_peri,boot_n);
    
    [opp_c,opp_tr]...
        = externalCBAnalysisGazeContingent_not_weighted_choice(weighted_signal_chosen,weighted_signal_notchosen,signal_chosen,signal_notchosen,params_boot(j,end-1),2,num_peri,boot_n);
    
    opp_cases(sub,:) = opp_c;
    opp_trials(sub) = opp_tr;
    sub_accuracy(sub) = accuracy;
    
    sub_opp_sacc(sub) = not_wt_trials_opp;
    trials_per_num_peri(sub) = (num_trials);
    num_frames = 1 + frames * 2;
    params_boot_fixed(sub,:,:) = params_boot;
    alpha(sub,:) = [prctile(params_boot(:, end).^2, 50) std(params_boot(:, end).^2)];
    bias(sub) = prctile(params_boot(:, end-1), 50);
    bias_err(sub,:) = [(prctile(params_boot(:, end-1), 50)-prctile(params_boot(:, end-1), 16)) (prctile(params_boot(:, end-1), 84)-prctile(params_boot(:, end-1), 50))];
    temporal_kernel(sub,:) = prctile(params_boot(:, 1:num_frames), 50);
    norm_temporal_kernel(sub,:) = temporal_kernel(sub,:)/mean(temporal_kernel(sub,:));
    lo_temporal_kernel(sub,:) = prctile(params_boot(:, 1:num_frames), 50) - prctile(params_boot(:, 1:num_frames), 16);
    hi_temporal_kernel(sub,:) = prctile(params_boot(:, 1:num_frames), 84) - prctile(params_boot(:, 1:num_frames), 50);
    %     beta(sub) = prctile(squeeze(abbl(:,2)),50);
    norm_all_linear(sub,:,:) = [sobl(:,1)/mean(temporal_kernel(sub,1:frames)) sobl(:,2)/mean(temporal_kernel(sub,1:frames)) sobl(:,3) sobl(:,4)];%sobl_norm;
    norm_slope_all(sub,:) = norm_all_linear(sub,:,1);
    norm_slope(sub) = prctile(squeeze(norm_all_linear(sub,:,1)),50);
    hprs_used(sub,:) = best_hprs;
    
    not_wt_belief_mid_opp(sub,:) = not_wt_belief_opp_cases(1,:);
    not_wt_lo_err_opp(sub,:) = not_wt_belief_opp_cases(3,:);
    not_wt_hi_err_opp(sub,:) = not_wt_belief_opp_cases(4,:);
    not_wt_prob_chose_in_favor_opp(sub,:) = not_wt_belief_opp_cases(2,:);
    
    
end


%%

txt = {'Ideal Observer (100 Samples)'; 'One Sample Approximation'; 'Two Samples Approximation'};
figure()
pr=2;
for i=1:(num_sub)
    disp(i)
    num_frames = frames;
    s(i)=subplot(2,num_sub,i);
    errorbar(1:num_frames,squeeze(temporal_kernel(i,1:frames)),squeeze(lo_temporal_kernel(i,1:frames)),squeeze(hi_temporal_kernel(i,1:frames)),'-ok','LineWidth',2)
    hold on;
    errorbar((1:num_frames-1),squeeze(temporal_kernel(i,(frames+1):(end - pr ))),squeeze(lo_temporal_kernel(i,(frames+1):(end - pr))),squeeze(hi_temporal_kernel(i,(frames+1):(end - pr))),'or','LineStyle','none')
    hold on;
    errorbar([num_frames num_frames]',squeeze(temporal_kernel(i,(end - pr+1 ):end)),squeeze(lo_temporal_kernel(i,(end - pr+1):end)),squeeze(hi_temporal_kernel(i,(end - pr+1):end)),'or','LineStyle','none')
    
    xlabel('Frames','Fontsize',15);
    ylabel('Norm Weights','Fontsize',15);
    title(txt{i});
    %     title(['Temporal PK for 2 in periphery for subject ' num2str(i)])
    hold('on');
    xlim([1 num_frames + 1]);
    text(1,0.8,['Trial Num: ' num2str(trials_per_num_peri(i))],'Fontsize',10);
    yline(0,'b');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    legend('chosen wts','unchosen wts','Location','southeast','Fontsize',10)
        ylim([-0.5 2.0]);
end


% figure(2)
for i=1:(num_sub)
    disp(i)
    s=subplot(2,num_sub,num_sub+i);
    errorbar(squeeze(not_wt_belief_mid_opp(i,:)),squeeze(not_wt_prob_chose_in_favor_opp(i,:)),squeeze(not_wt_lo_err_opp(i,:)),squeeze(not_wt_hi_err_opp(i,:)),'-bo','LineWidth',2);
    hold on;
    yline(0.5,'k','LineWidth',2);
    xlabel('Accumulated evidence','Fontsize',15)
    ylabel('Prob chose in favor','Fontsize',15)
    title(txt{i});
    %     title(['Ext CB for 2 in periphery for subject ' num2str(i)])
        xlim([0 500]);
        ylim([0 1.5]);
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    text(s.XLim(1),s.YLim(2)-0.15,['Opp Saccade Num: ' num2str(floor(mean(sub_opp_sacc(i))))],'Fontsize',10);
end


%%

figure(6)
subplot(1,3,1)
hold on;
bar(1:num_sub,sub_accuracy./trials);
hold on;
% xticks(txt);
yline(0.5,'k');
xlabel('Subject Number','Fontsize',20)
ylabel('Percent Correct','Fontsize',20)
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
% title(['Accuracy under fixed noise and ratio of 0.7'],'Fontsize',18)


% this figure plots the accuracy of the subjects
% figure(7)
subplot(1,3,2)
hold on;
bar(1,opp_cases(2,2),'g');
hold on;
bar(2,opp_cases(3,2),'r');
hold on;
bar(3,opp_cases(1,2),'b');
% bar(1:num_sub,squeeze(opp_cases([2 3 1],2)),);
hold on;
% xticks(txt);
errorbar(1:num_sub,squeeze(opp_cases([2 3 1],2)),squeeze(opp_cases([2 3 1],3)),squeeze(opp_cases([2 3 1],4)),'ok')
xlabel('Subject Number','Fontsize',20)
ylabel('Saccade probability to evidence in favor of already accumulated evidence','Fontsize',20)
yline(0.5,'k');
ylim([0.4 1]);
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
% xlim([1 3]);
xticklabels({'B','MB','I'});
% title(['Collapsed Prob chose in favor of accumulated evidence'],'Fontsize',18)

%%
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

