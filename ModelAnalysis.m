clear all; close all;clc;
simulated_data = load('summary_simulation_3.mat');
num_sub = length(simulated_data.subj);%testing simulated subjects with various degrees of approximation
num_peri = 2;%stimuli in periphery is 2
frames = 4;% number of fixations, saccades=frames-1
num_images = num_peri * frames + 1;%total number of stimuli present on the screen during a trial
opt_type = 'model'; % type pf optimization in PK computation
hpr1 = 0.0;
hpr2 = 10;%logspace(-1, 5, 7);
points = [1 2 5];
boot_n = 5;

for pp=1:length(points)
    not_wt_belief_extCB_opp = zeros(num_sub,points(pp));
    not_wt_belief_err_opp = zeros(num_sub,points(pp));
    not_wt_lobelief_err_opp = zeros(num_sub,points(pp));
    not_wt_hibelief_err_opp = zeros(num_sub,points(pp));
    not_wt_probability_chose_in_favor_opp = zeros(num_sub,points(pp));
end

for sub=1:num_sub
    disp(['Running analysis for simulated subject ' num2str(sub) ' ...'])
    signal_all = simulated_data.subj{sub}.signals;
    trials = size(simulated_data.subj{sub}.signals,1);
    choice_raw = simulated_data.subj{sub}.choice;
    choice_raw(choice_raw==-1) = 0;
    disp(['Data obtained for subject ' num2str(sub)]);
    %     % this part computes the hyperparameters for the PK
    [best_hprs, ~] = CustomRegression.xValidatePK_with_lapse(signal_all, choice_raw, frames, hpr1, 0, hpr2, 0, 10, opt_type);
    disp('Best hyperparameters found!');
    % this part computes the PK weights for the chosen and the not chosen
    % signals, to be later used for computation of accumulated evidence
    all_frames = 1 + frames * num_peri;
    for j = 1:boot_n
        [signal, choice] = bootstrap(signal_all, choice_raw, trials);
        [sobl(j,:), ~] =  CustomRegression.LinearPK_with_lapse(signal(:,1:frames), choice, 0);
        if mod(j,100)==0 || j==1
            disp(['Bootstrap number ' num2str(j) ' ...']);
        end
        [params_boot(j,:), ~, ~, ~, ~, ~] =  CustomRegression.PsychophysicalKernelwithlapse(signal, choice, frames, best_hprs(1), 0, best_hprs(3), 0, opt_type);
        [abbl(j,:), ~, ~] =  CustomRegression.ExponentialPK_with_lapse(signal(:,1:frames), choice, 0);
    end
    
    % signals weighted by PK weights based on final categorical choice
    weighted_signal_all = prctile(params_boot(:, 1:all_frames), 50).* signal_all;
    
    weighted_signal_chosen = weighted_signal_all(:,1:frames);
    weighted_signal_notchosen = weighted_signal_all(:,frames+1:end-num_peri);
    weighted_signal_notchosen = reshape(weighted_signal_notchosen,trials,num_peri-1,frames-1);
    
    % ideal signals not weighted at all
    signal_chosen = signal_all(:,1:frames);
    signal_notchosen = signal_all(:,frames+1:end-num_peri);
    signal_notchosen = reshape(signal_notchosen,trials,num_peri-1,frames-1);
    
    disp('Computing bias in saccade ...');
    % computing the effects of external CB
    for pp=1:length(points)
        [not_wt_belief_opp_cases{pp},not_wt_trials_opp{pp}]...
            = externalCBAnalysisGazeContingent_not_weighted_choice(weighted_signal_chosen,weighted_signal_notchosen,signal_chosen,signal_notchosen,params_boot(j,end-1),points(pp)+1,num_peri,boot_n);
        not_wt_belief_mid_opp{pp}(sub,:) = not_wt_belief_opp_cases{pp}(1,:);
        not_wt_lo_err_opp{pp}(sub,:) = not_wt_belief_opp_cases{pp}(3,:);
        not_wt_hi_err_opp{pp}(sub,:) = not_wt_belief_opp_cases{pp}(4,:);
        not_wt_prob_chose_in_favor_opp{pp}(sub,:) = not_wt_belief_opp_cases{pp}(2,:);
    end
    
    sub_accuracy(sub) = sum(simulated_data.subj{1}.choice==simulated_data.subj{1}.correct)/trials;
    sub_opp_sacc(sub) = not_wt_trials_opp{1};
    trials_per_num_peri(sub) = (trials);
    num_frames = all_frames;
    
    params_boot_fixed(sub,:,:) = params_boot;
    alpha(sub,:) = [prctile(params_boot(:, end).^2, 50) std(params_boot(:, end).^2)];
    bias(sub) = prctile(params_boot(:, end-1), 50);
    bias_err(sub,:) = [(prctile(params_boot(:, end-1), 50)-prctile(params_boot(:, end-1), 16)) (prctile(params_boot(:, end-1), 84)-prctile(params_boot(:, end-1), 50))];
    temporal_kernel(sub,:) = prctile(params_boot(:, 1:num_frames), 50);
    norm_temporal_kernel(sub,:) = temporal_kernel(sub,:)/mean(temporal_kernel(sub,:));
    lo_temporal_kernel(sub,:) = prctile(params_boot(:, 1:num_frames), 50) - prctile(params_boot(:, 1:num_frames), 16);
    hi_temporal_kernel(sub,:) = prctile(params_boot(:, 1:num_frames), 84) - prctile(params_boot(:, 1:num_frames), 50);
    beta(sub) = prctile(squeeze(abbl(:,2)),50);
    norm_all_linear(sub,:,:) = [sobl(:,1)/mean(temporal_kernel(sub,1:frames)) sobl(:,2)/mean(temporal_kernel(sub,1:frames)) sobl(:,3) sobl(:,4)];%sobl_norm;
    norm_slope_all(sub,:) = norm_all_linear(sub,:,1);
    norm_slope(sub) = prctile(squeeze(norm_all_linear(sub,:,1)),50);
    hprs_used(sub,:) = best_hprs;
    
    disp(['Analysis for simulated subject ' num2str(sub) ' complete!!']);
    disp('-----------------------------------------------------------------------------------------------------');
end


%%
for i=1:num_sub
    if simulated_data.subj{i}.nsamp==100
        txt{i} = ['Ideal(' num2str(simulated_data.subj{i}.nsamp) ' samples)'];
    elseif simulated_data.subj{i}.nsamp==1
        txt{i} = [num2str(simulated_data.subj{i}.nsamp) ' sample'];
    else
        txt{i} = [num2str(simulated_data.subj{i}.nsamp) ' samples'];
    end
end
figure();
pr = 2;
for i=1:(num_sub)
    num_frames = frames;
    s(i) = subplot(2,num_sub,i);
    errorbar(1:num_frames,squeeze(temporal_kernel(i,1:frames)),squeeze(lo_temporal_kernel(i,1:frames)),squeeze(hi_temporal_kernel(i,1:frames)),'-ok','LineWidth',2);
    hold on;
    errorbar((2:num_frames),squeeze(temporal_kernel(i,(frames+1):(end - pr ))),squeeze(lo_temporal_kernel(i,(frames+1):(end - pr))),squeeze(hi_temporal_kernel(i,(frames+1):(end - pr))),'or','LineStyle','none');
    hold on;
    errorbar([num_frames+1 num_frames+1]',squeeze(temporal_kernel(i,(end - pr+1 ):end)),squeeze(lo_temporal_kernel(i,(end - pr+1):end)),squeeze(hi_temporal_kernel(i,(end - pr+1):end)),'or','LineStyle','none');
    if i==2
        xlabel('Fixations','Fontsize',20);
    end
    if i==1
        ylabel('Weights','Fontsize',20);
    end
    xticks([1 2 3 4]);
    title(txt{i});
    hold('on');
    xlim([1 num_frames + 1]);
    text(1,0.8,['Trial Num: ' num2str(trials_per_num_peri(i))],'Fontsize',10);
    yline(0,'b');
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    legend('Chosen wts','Unchosen wts','Location','southeast','Fontsize',10)
    ylim([-0.1 0.1]);
end

for i=1:(num_sub)
    s = subplot(2,num_sub,num_sub+i);
    errorbar(squeeze(not_wt_belief_mid_opp{3}(i,:)),squeeze(not_wt_prob_chose_in_favor_opp{3}(i,:)),squeeze(not_wt_lo_err_opp{3}(i,:)),squeeze(not_wt_hi_err_opp{3}(i,:)),'-bo','LineWidth',2);
    yline(0.5,'k','LineWidth',2);
    if i==2
        xlabel('Accumulated evidence','Fontsize',15);
    end
    if i==1
        ylabel('Prob. chose in favor','Fontsize',15);
    end
    title(txt{i});
    xlim([0 50]);
    ylim([0.4 1.0]);
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    text(s.XLim(1),s.YLim(2)-0.15,['Opp Saccade Num: ' num2str(floor(mean(sub_opp_sacc(i))))],'Fontsize',10);
end
%%
for i=1:num_sub
    if simulated_data.subj{i}.nsamp==1
        txt{i} = [num2str(simulated_data.subj{i}.nsamp) ' sample'];
    else
        txt{i} = [num2str(simulated_data.subj{i}.nsamp) ' samples'];
    end
end
figure();
clr = ['g', 'r', 'b'];
model_type = {'B', 'MB', 'I'};
pr = 2;
subplot(2,5,5)
for i=1:(num_sub)
    num_frames = frames;
    LH(i) = errorbar(1:num_frames,squeeze(temporal_kernel(i,1:frames)),squeeze(lo_temporal_kernel(i,1:frames)),squeeze(hi_temporal_kernel(i,1:frames)),['-o' clr(i)],'LineWidth',2);
    L{i} = txt{i};%model_type{i};
    hold on;
    errorbar((2:num_frames),squeeze(temporal_kernel(i,(frames+1):(end - pr ))),squeeze(lo_temporal_kernel(i,(frames+1):(end - pr))),squeeze(hi_temporal_kernel(i,(frames+1):(end - pr))),'Marker','o','Markersize',8,'MarkerFaceColor',clr(i),'MarkerEdgeColor','k','Linestyle','none');
    hold on;
    errorbar([num_frames+1 num_frames+1]',squeeze(temporal_kernel(i,(end - pr+1 ):end)),squeeze(lo_temporal_kernel(i,(end - pr+1):end)),squeeze(hi_temporal_kernel(i,(end - pr+1):end)),'Marker','o','Markersize',8,'MarkerFaceColor',clr(i),'MarkerEdgeColor','k','Linestyle','none');
    xlabel('Fixations','Fontsize',15);
    ylabel('Weights','Fontsize',15);
    hold('on');
    xlim([1 num_frames + 1]);
    yline(0,'k','linewidth',2);
    ax = gca;
    set(ax, 'box','off');
    ax.LineWidth=2;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    ylim([-0.01 0.1]);
    xticks([1 2 3 4]);
end
legend(LH,L,'Fontsize',10, 'Box','off');

subplot(2,5,3)
for i=1:(num_sub)
    errorbar(squeeze(not_wt_belief_mid_opp{2}(i,:)),squeeze(not_wt_prob_chose_in_favor_opp{2}(i,:)),squeeze(not_wt_lo_err_opp{2}(i,:)),squeeze(not_wt_hi_err_opp{2}(i,:)),['-o' clr(i)],'LineWidth',2);
    hold on;
    yline(0.5,'k','LineWidth',2);
    xlabel({'Accumulated';' evidence'},'Fontsize',15);
    ylabel('Prob chose in favor','Fontsize',15);
    xlim([0 30]);
    xticks([0 15 30]);
    ylim([0.45 1.0]);
    ax = gca;
    set(ax, 'box','off');
    ax.LineWidth=2;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
end

subplot(2,5,4)
for i=1:(num_sub)
    errorbar(squeeze(not_wt_belief_mid_opp{3}(i,:)),squeeze(not_wt_prob_chose_in_favor_opp{3}(i,:)),squeeze(not_wt_lo_err_opp{3}(i,:)),squeeze(not_wt_hi_err_opp{3}(i,:)),['-o' clr(i)],'LineWidth',2);
    hold on;
    yline(0.5,'k','LineWidth',2);
    xlabel({'Accumulated';' evidence'},'Fontsize',15);
    ylabel('Prob chose in favor','Fontsize',15);
    xlim([0 30]);
    xticks([0 15 30]);
    ylim([0.4 1.0]);
    ax = gca;
    set(ax, 'box','off');
    ax.LineWidth=2;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
end

subplot(2,5,1)
for mb=1:length(model_type)
    bar(mb,sub_accuracy(mb),clr(mb));
    hold on;
end
hold on;
errorbar(1:num_sub,sub_accuracy,sqrt((sub_accuracy.*(1-sub_accuracy))./trials),'ok','linewidth',2);
hold on;
yline(0.5,'k','linewidth',2);
ylim([0.4 1.0]);
xlabel('Subject Type','Fontsize',20);
ylabel('Percent Correct','Fontsize',20);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
xlim([0.25 3.75]);
xticks(linspace(1,length(model_type),length(model_type)));
xticklabels({'B','MB','I'});

subplot(2,5,2)
for mb=1:length(model_type)
    bar(mb,squeeze(not_wt_prob_chose_in_favor_opp{1}(mb)),clr(mb));
    hold on;
end
hold on;
errorbar(1:num_sub,squeeze(not_wt_prob_chose_in_favor_opp{1}(:)),squeeze(not_wt_lo_err_opp{1}(:)),squeeze(not_wt_hi_err_opp{1}(:)),'ok','linewidth',2);
xlabel('Subject Type','Fontsize',20)
ylabel({'Prob chose in favor'},'Fontsize',20)
yline(0.5,'k','linewidth',2);
ylim([0.45 1]);
ax = gca;
set(ax, 'box','off');
ax.LineWidth = 2;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
xlim([0.25 3.75]);
xticks(linspace(1,length(model_type),length(model_type)));
xticklabels({'B','MB','I'});

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

