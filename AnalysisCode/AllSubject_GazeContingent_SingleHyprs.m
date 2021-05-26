% This is a script to display results of the analysis for the external CB project
% initialize parameters and hyperparameters
clear all; close all; clc;
warning('off');
sanity_bins = 10;
boots = 500;
points = [1 2 3]; % number of points in accumulated belief vs bias plot
ps_points = 3;
hpr_ridge = logspace(-3, 3, 14);
hpr_ar1 = 0;
hpr_curvature = logspace(-3, 3, 14);
standardize = 0;
folds = 10;
frames = 4; % number of fixations
expt_type = 3;

subjects = ...
    {
    'GCV3-subject01'; 'GCV3-subject02';...
    'GCV3-subject03'; 'GCV3-subject07';...
    'GCV3-subject08'; 'GCV3-subject11';...
    'GCV3-subject13'; 'GCV3-subject15';...
    'GCV3-subject16'; 'GCV3-subject17';...
    'GCV3-subject18'; 'GCV3-subject19';...
    'GCV3-subject20'; 'GCV3-subject21';...
    'GCV3-subject22'; 'GCV3-subject23';...
    };


disp('Starting search of best hyper parameters across all subjects...');
[num_sub,~] = size(subjects);
% [best_hprs] = CustomRegression.combined_cross_validation_hprs(subjects, expt_type, hpr_ridge, hpr_ar1, hpr_curvature, standardize, folds);
best_hprs = [0.2955         0    0.4437];

fx = @(mu,nt) normcdf(abs(mu-0.5)./sqrt((mu*(1-mu))/nt),'upper');

% pass on the parameters and subject ids to do respective analysis for each case of num_peri(number of elements in the periphery of the experiment: 2)
for sub=1:num_sub
    disp(['Starting analysis for subject ' num2str(sub) '/' num2str(num_sub) ' ...']);
    tic;
    [params_boot, sobl, abbl, log_bernoulli{sub},...
        data, accuracy, belief_opp_cases, ...
        not_wt_trials_opp,opp_r,...
        mn_sig(sub,:),mn_sig_norm(sub,:),mn_sig_ideal_foveal_only(sub,:),mn_sig_ideal_foveal_only_norm(sub,:),...
        mn_sig_ideal_all(sub,:),mn_sig_ideal_all_norm(sub,:),mn_perf(sub,:),mn_perf_ideal_foveal_only(sub,:),mn_perf_ideal_all(sub,:),...
        tr_seg(sub,:),err_perf(sub,:),tr_seg_ideal_foveal_only(sub,:),err_perf_ideal_foveal_only(sub,:),tr_seg_ideal_all(sub,:),err_perf_ideal_all(sub,:),...
        mn_sig_ideal_all_actual_sig(sub,:),mn_perf_actual_sig(sub,:),mn_perf_ideal_foveal_only_actual_sig(sub,:),mn_perf_ideal_all_actual_sig(sub,:),...
        tr_seg_actual_sig(sub,:),err_perf_actual_sig(sub,:),tr_seg_ideal_foveal_only_actual_sig(sub,:),err_perf_ideal_foveal_only_actual_sig(sub,:),...
        tr_seg_ideal_all_actual_sig(sub,:),err_perf_ideal_all_actual_sig(sub,:),...
        sanity_actual_data(sub,:,:), sanity_predicted_data(sub,:,:), sanity_mean_sig_shown(sub,:), prob_opp_divided(sub,:,:),...
        belief_opp_cases_first_half,belief_opp_cases_second_half,accuracy_first_half,accuracy_second_half]...
        = AnalysisGazeContingent_SingleHyprs(subjects{sub}, expt_type, boots, points+1, best_hprs, ps_points, sanity_bins, standardize);
    
    choice{sub} = data.choice;
    %assign PK weight related values from the returned results
    sub_opp_sacc(sub) = not_wt_trials_opp{1};
    trials_per_num_peri(sub) = (sum(data.num_peri(1:data.current_trial)==2));
    num_frames = 1 + frames * 2;
    %     alpha(sub,:) = [prctile(params_boot(:,end).^2,50) sqrt(var(params_boot(:,end).^2)/size(params_boot,1))];%[prctile(params_boot(:, end).^2, 50) std(params_boot(:, end).^2)];
    alpha(sub,:) = [prctile(1e-4+(1-1e-4)*sigmoid(params_boot(:,end)),50) sqrt(var(1e-4+(1-1e-4)*sigmoid(params_boot(:,end)))/size(params_boot,1))];%[prctile(params_boot(:, end).^2, 50) std(params_boot(:, end).^2)];
    bias(sub) = prctile(params_boot(:, end-1), 50);
    temporal_kernel(sub,:) = prctile(params_boot(:, 1:num_frames), 50);
    norm_temporal_kernel(sub,:) = temporal_kernel(sub,:)/mean(temporal_kernel(sub,:));
    lo_temporal_kernel(sub,:) = prctile(params_boot(:, 1:num_frames), 50) - prctile(params_boot(:, 1:num_frames), 16);
    hi_temporal_kernel(sub,:) = prctile(params_boot(:, 1:num_frames), 84) - prctile(params_boot(:, 1:num_frames), 50);
    beta(sub) = prctile(squeeze(abbl(:,2)),50);
    norm_all_linear(sub,:,:) = [sobl(:,1)/mean(temporal_kernel(sub,1:frames)) sobl(:,2)/mean(temporal_kernel(sub,1:frames)) sobl(:,3) sobl(:,4)];%sobl_norm;
    norm_slope_all(sub,:) = norm_all_linear(sub,:,1);
    norm_slope(sub) = prctile(squeeze(norm_all_linear(sub,:,1)),50);
    hprs_used(sub,:) = best_hprs;
    
    for pp=1:length(points)
        belief_mid_opp{pp}(sub,:) = belief_opp_cases{pp}(1,:);
        lo_err_opp{pp}(sub,:) = belief_opp_cases{pp}(3,:);
        hi_err_opp{pp}(sub,:) = belief_opp_cases{pp}(4,:);
        prob_chose_in_favor_opp{pp}(sub,:) = belief_opp_cases{pp}(2,:);
        opp_raw_sig{pp}(sub,:) = opp_r{pp}(1,:);
        opp_raw_prob{pp}(sub,:) = opp_r{pp}(2,:);
        opp_raw_lo_err{pp}(sub,:) = opp_r{pp}(3,:);
        opp_raw_hi_err{pp}(sub,:) = opp_r{pp}(4,:);
    end
    belief_mid_opp_first_half(sub,:) = belief_opp_cases_first_half(1,:);
    lo_err_opp_first_half(sub,:) = belief_opp_cases_first_half(3,:);
    hi_err_opp_first_half(sub,:) = belief_opp_cases_first_half(4,:);
    prob_chose_in_favor_opp_first_half(sub,:) = belief_opp_cases_first_half(2,:);
    belief_mid_opp_second_half(sub,:) = belief_opp_cases_second_half(1,:);
    lo_err_opp_second_half(sub,:) = belief_opp_cases_second_half(3,:);
    hi_err_opp_second_half(sub,:) = belief_opp_cases_second_half(4,:);
    prob_chose_in_favor_opp_second_half(sub,:) = belief_opp_cases_second_half(2,:);
    
    opp_trials(sub) = not_wt_trials_opp{1};
    sub_accuracy(sub) = sum(accuracy);
    sub_accuracy_first_half(sub) = accuracy_first_half;
    sub_accuracy_second_half(sub) = accuracy_second_half;
    
    p_values_opp_cases(sub) = fx(prob_chose_in_favor_opp{1}(sub,:),sub_opp_sacc(sub));
    if (p_values_opp_cases(sub))<0.001
        stars{sub} = '***';
    elseif (p_values_opp_cases(sub))<0.01
        stars{sub} = '**';
    elseif (p_values_opp_cases(sub))<0.05
        stars{sub} = '*';
    else
        stars{sub} = ' ';
    end
    disp(['The learned regression weights: ' num2str(temporal_kernel(sub,:))]);
    disp(['The learned lapse: ' num2str(alpha(sub,1))]);
    disp(['The learned bias: ' num2str(bias(sub))]);
    disp(['The significance level: ' stars{sub}]);
    disp(['Starting pattern matching analysis for subject ' num2str(sub) ' ...']);
    
    [pattern_one_previous{sub},belief_one_previous{sub},...
        pattern_two_previous{sub}, belief_two_previous{sub},...
        pattern_three_previous{sub}, belief_three_previous{sub}]...
        = pattern_matching_comparison(data,squeeze(temporal_kernel(sub,:)),bias(sub));
    
    [pattern_one_previous_case_belief_evi_opp{sub},belief_one_previous_case_belief_evi_opp{sub}]...
    = pattern_matching_comparison_belief_opp_fovea(data,squeeze(temporal_kernel(sub,:)),bias(sub));
    
    disp(['Completed analysis for subject ' num2str(sub) ' !!']);
    toc;
    disp('-----------------------------------------------------------------------------------------------------');
    
    figure(sub);
    subplot(1,2,1)
    pr = 2;
    errorbar(1:frames,squeeze(temporal_kernel(sub,1:frames)),squeeze(lo_temporal_kernel(sub,1:frames)),squeeze(hi_temporal_kernel(sub,1:frames)),'-ob','LineWidth',2);
    hold on;
    errorbar((2:frames),squeeze(temporal_kernel(sub,(frames+1):(end - pr ))),squeeze(lo_temporal_kernel(sub,(frames+1):(end - pr))),squeeze(hi_temporal_kernel(sub,(frames+1):(end - pr))),'or','LineStyle','none','LineWidth',2);
    hold on;
    errorbar([frames+1 frames+1]',squeeze(temporal_kernel(sub,(end - pr+1 ):end)),squeeze(lo_temporal_kernel(sub,(end - pr+1):end)),squeeze(hi_temporal_kernel(sub,(end - pr+1):end)),'or','LineStyle','none','LineWidth',2);
    xlabel('Fixation Number','fontsize',30);
    ylabel('Unnormalized Weights','fontsize',30);
    hold('on');
    xlim([1 frames + 1]);
    yline(0,'k','linewidth',2);
    %     ylim([-0.2 0.35]);
    xticks([1 2 3 4]);
    ax = gca;
    set(ax, 'box','off');
    ax.LineWidth = 2;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    hold on;
    
    subplot(1,2,2)
    errorbar(squeeze(sanity_mean_sig_shown(sub,:)),squeeze(sanity_actual_data(sub,:,1)),squeeze(sanity_actual_data(sub,:,2)),'-ob','LineWidth',2);
    hold on;
    errorbar(squeeze(sanity_mean_sig_shown(sub,:)),squeeze(sanity_predicted_data(sub,:,1)),squeeze(sanity_predicted_data(sub,:,2)),'-or','LineWidth',2);
    xlabel('Signal Displayed','fontsize',30);
    ylabel('Probability chose left','fontsize',30);
    hold('on');
    yline(0.5,'k','linewidth',2);
    xline(0.0,'k','linewidth',2);
    ylim([0.0 1.0]);
    ax = gca;
    set(ax, 'box','off');
    ax.LineWidth = 2;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    %     pause;
    sgtitle(['Analysis for Subject ' num2str(sub)],'fontsize',30);
end

if mod(num_sub,2)==0
    sb_plt = num_sub/2;
    md_x = floor(sb_plt/2);
else
    sb_plt = floor(num_sub/2) + 1;
    md_x = floor(sb_plt/2) + 1;
end
authors = [2 3];
chosen_sub = [];
for i=1:num_sub
    if stars{i}==' '
        chosen_sub(end+1) = i;
    end
end
    
%%
fig2 = figure();
set(fig2,'defaultLegendAutoUpdate','off');
for sub=1:(num_sub)
    s(sub)=subplot(2,sb_plt,sub);
    if sub==1 || sub==(sb_plt+1)
        LH(1) = errorbar(squeeze(sanity_mean_sig_shown(sub,:)),squeeze(sanity_actual_data(sub,:,1)),squeeze(sanity_actual_data(sub,:,2)),'-ob','LineWidth',2);
        L{1} = 'Data';
        hold on;
        LH(2) = errorbar(squeeze(sanity_mean_sig_shown(sub,:)),squeeze(sanity_predicted_data(sub,:,1)),squeeze(sanity_predicted_data(sub,:,2)),'-or','LineWidth',2);
        L{2} = 'Predicted';
        hold on;
    else
        errorbar(squeeze(sanity_mean_sig_shown(sub,:)),squeeze(sanity_actual_data(sub,:,1)),squeeze(sanity_actual_data(sub,:,2)),'-ob','LineWidth',2);
        hold on;
        errorbar(squeeze(sanity_mean_sig_shown(sub,:)),squeeze(sanity_predicted_data(sub,:,1)),squeeze(sanity_predicted_data(sub,:,2)),'-or','LineWidth',2);
        hold on;
    end
    if sub==md_x || sub==(sb_plt+md_x)
        xlabel('Signal Displayed','fontsize',30);
    end
    if sub==1 || sub==(sb_plt+1)
        ylabel('Probability chose left','fontsize',30);
    end
    if sub==1
        legend(LH,L, 'Fontsize',20, 'Box','off','location','northwest');
    end
    hold('on');
    text(200,0.25,[num2str(trials_per_num_peri(sub)) ' trials'],'fontsize',15);
    yline(0.5,'k','linewidth',2);
    xline(0,'k','linewidth',2);
    ylim([0.0 1.0]);
    ax = gca;
    set(ax, 'box','off');
    ax.LineWidth=2;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title(['Subject ' num2str(sub)],'fontsize',20);
end
sgtitle('Checking prediction w.r.t fits','fontsize',30);

%%
fig = figure();
set(fig,'defaultLegendAutoUpdate','off');
pr=2;
for i=1:(num_sub)
    s(i)=subplot(2,sb_plt,i);
    if i==1 || i==(sb_plt+1)
        LH(1)=errorbar(1:frames,squeeze(temporal_kernel(i,1:frames)),squeeze(lo_temporal_kernel(i,1:frames)),squeeze(hi_temporal_kernel(i,1:frames)),'-ob','LineWidth',2);
        L{1} = 'Foveated';
        hold on;
        LH(2)= errorbar((2:frames),squeeze(temporal_kernel(i,(frames+1):(end - pr ))),squeeze(lo_temporal_kernel(i,(frames+1):(end - pr))),squeeze(hi_temporal_kernel(i,(frames+1):(end - pr))),'or','LineStyle','none','LineWidth',2);
        L{2} = 'Not foveated';
        hold on;
    else
        errorbar(1:frames,squeeze(temporal_kernel(i,1:frames)),squeeze(lo_temporal_kernel(i,1:frames)),squeeze(hi_temporal_kernel(i,1:frames)),'-ob','LineWidth',2);
        hold on;
        errorbar((1:frames-1),squeeze(temporal_kernel(i,(frames+1):(end - pr ))),squeeze(lo_temporal_kernel(i,(frames+1):(end - pr))),squeeze(hi_temporal_kernel(i,(frames+1):(end - pr))),'or','LineStyle','none','LineWidth',2);
        hold on;
        errorbar([frames frames]',squeeze(temporal_kernel(i,(end - pr+1 ):end)),squeeze(lo_temporal_kernel(i,(end - pr+1):end)),squeeze(hi_temporal_kernel(i,(end - pr+1):end)),'or','LineStyle','none','LineWidth',2);
        hold on;
    end
    if i==md_x || i==(sb_plt+md_x)
        xlabel('Fixation Number','fontsize',30);
        hold on;
    end
    if i==1 || i==(sb_plt+1)
        ylabel('Unnormalized Weights','fontsize',30);
        hold on;
    end
    if i==1
        legend(LH,L, 'Fontsize',20, 'Box','off');
    end
    hold on;
    xlim([1 frames + 1]);
    text(2,2,[num2str(trials_per_num_peri(i)) ' trials'],'fontsize',15);
    yline(0,'k','linewidth',2);
    ylim([-2 5]);
    xticks([1 2 3 4]);
    ax = gca;
    set(ax, 'box','off');
    ax.LineWidth=2;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title(['Subject ' num2str(i)],'fontsize',20);
end
sgtitle('Learned regression weights','fontsize',30);
%%
figure();
for i=1:(num_sub)
    s=subplot(2,sb_plt,i);
    bar(1,squeeze(prob_chose_in_favor_opp{1}(i,:)),'FaceColor','b','EdgeColor','k','LineWidth',0.75);
    hold on;
    errorbar(1,squeeze(prob_chose_in_favor_opp{1}(i,:)),squeeze(lo_err_opp{1}(i,:)),squeeze(hi_err_opp{1}(i,:)),'-ko','LineWidth',2,'linestyle','none');
    hold on;
    yline(0.5,'k','LineWidth',2);
    if i==1 || i==(sb_plt+1)
        ylabel('Probability','fontsize',30)
    end
    perf = sub_accuracy(i)/trials_per_num_peri(i);
    bar(2,perf,'FaceColor','r','EdgeColor','k','LineWidth',0.75);
    hold on;
    errorbar(2,sub_accuracy(i)/trials_per_num_peri(i),sqrt((perf*(1-perf))/trials_per_num_peri(i)),'-ko','LineWidth',2,'linestyle','none');
    hold on;
    yline(0.5,'k','LineWidth',2);
    xlim([0.25 2.75]);
    ylim([0.3 1.0]);
    text(0.5,0.9,['Opp. Sacc: ' num2str(floor(mean(sub_opp_sacc(i))))],'fontsize',15);
    xticks([1 2]);
    xticklabels({'Bias', 'Perf.'});
    ax = gca;
    set(ax, 'box','off');
    ax.LineWidth=2;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title(['Subject ' num2str(i)],'fontsize',20);
end
sgtitle('Subjectwise performance and bias','fontsize',30);
%%
for pp=2:length(points)
    figure();
    for i=1:(num_sub)
        s=subplot(2,sb_plt,i);
        errorbar(squeeze(belief_mid_opp{pp}(i,:)),squeeze(prob_chose_in_favor_opp{pp}(i,:)),squeeze(lo_err_opp{pp}(i,:)),squeeze(hi_err_opp{pp}(i,:)),'-bo','LineWidth',2);
        hold on;
        %     errorbar(squeeze(opp_raw_sig{pp}(i,:)),squeeze(opp_raw_prob{pp}(i,:)),squeeze(opp_raw_lo_err{pp}(i,:)),squeeze(opp_raw_hi_err{pp}(i,:)),'-ro','LineWidth',2);
        yline(0.5,'k','LineWidth',2);
        if i==1 || i==(sb_plt+1)
            ylabel('Probability chose in favor','fontsize',30);
        end
        if i==md_x || i==(sb_plt+md_x)
            xlabel('Accumulated evidence','fontsize',30);
        end
        yline(0.5,'k','LineWidth',2);
        xticks(round(squeeze(belief_mid_opp{pp}(i,:)),2));
        xlim([min(belief_mid_opp{pp}(i,:))-0.01 max(squeeze(belief_mid_opp{pp}(i,:)))+0.01]);
        ylim([0.0 1.0]);
        text(belief_mid_opp{pp}(i,1)+0.1,0.2,['Opp. Sacc: ' num2str(floor(mean(sub_opp_sacc(i))))],'fontsize',15);
        ax = gca;
        set(ax, 'box','off');
        ax.LineWidth=2;
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        title(['Subject ' num2str(i)],'fontsize',20);
    end
    sgtitle(['Bias w.r.t accumulated evidence for ' num2str(pp) ' bins'],'fontsize',30);
end

%%
figure();
for i=1:(num_sub)
    s=subplot(2,sb_plt,i);
    errorbar([1 2 3],squeeze(prob_opp_divided(i,:,1)),squeeze(prob_opp_divided(i,:,2)),'-bo','LineWidth',2);
    hold on;
    yline(0.5,'k','LineWidth',2);
    if i==1 || i==(sb_plt+1)
        ylabel('Probability chose in favor','fontsize',30);
    end
    if i==md_x || i==(sb_plt+md_x)
        xlabel('Saccade Number or Time','fontsize',30);
    end
    yline(0.5,'k','LineWidth',2);
    xticks([1 2 3]);
    xlim([0.75 3.25]);
    ylim([0.2 1.0]);
    ax = gca;
    set(ax, 'box','off');
    ax.LineWidth=2;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title(['Subject ' num2str(i)],'fontsize',20);
end
sgtitle(['Bias w.r.t Time or Saccade Number'],'fontsize',30);

figure();
for i=1:(num_sub)
    errorbar([1 2 3],squeeze(prob_opp_divided(i,:,1)),2*squeeze(prob_opp_divided(i,:,2)),'--bo','LineWidth',1);
    hold on;
    yline(0.5,'k','LineWidth',2);
    ylabel('Probability chose in favor','fontsize',30);
    xlabel('Saccade Number or Time','fontsize',30);
   
    yline(0.5,'k','LineWidth',2);
    xticks([1 2 3]);
    xlim([0.75 3.25]);
    ylim([0.2 1.0]);
    ax = gca;
    set(ax, 'box','off');
    ax.LineWidth=2;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
end
plot([1 2 3],mean(squeeze(prob_opp_divided(:,:,1)),1),'-ko','LineWidth',5);
sgtitle(['Bias w.r.t Time or Saccade Number'],'fontsize',30);


%%
figure();
subplot(2,3,1)
hold on;
bar(1:num_sub,sub_accuracy./trials_per_num_peri, 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',0.75);
hold on;
bar(authors,sub_accuracy(authors)./trials_per_num_peri(authors), 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor',[0.6350 0.0780 0.1840],'LineWidth',3.0);
hold on;
yline(0.5,'k','LineWidth',1.5);
xlabel('Subject Number','Fontsize',20);
ylabel('Percent Correct','Fontsize',20);
ax = gca;
ax.LineWidth=2;
ylim([0.4 1]);
xlim([0.5 num_sub+0.5]);
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
ax.XTick = linspace(1,num_sub,num_sub);
hold on;
tmp = sub_accuracy./trials_per_num_peri;
errorbar(1:num_sub,tmp,sqrt((tmp.*(1-tmp))./trials_per_num_peri),'ok','LineWidth', 1.5);
hold on;

subplot(2,3,2)
hold on;
bar(1:num_sub,squeeze(prob_chose_in_favor_opp{1}), 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',0.75);
bar(authors,squeeze(prob_chose_in_favor_opp{1}(authors)), 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor',[0.6350 0.0780 0.1840],'LineWidth',3);
hold on;
errorbar(1:num_sub,squeeze(prob_chose_in_favor_opp{1}),squeeze(lo_err_opp{1}),squeeze(hi_err_opp{1}),'ok','LineWidth', 1.5)
hold on;
text([1:num_sub]-0.3,squeeze(prob_chose_in_favor_opp{1})+0.065,stars,'FontSize',15,'FontWeight','bold');
xlabel('Subject Number','Fontsize',20)
ylabel({'Probability of saccading to ';'confirming evidence'},'Fontsize',20);
yline(0.5,'k','LineWidth',1.5);
ax = gca;
ax.LineWidth=2;
ylim([0.4 1]);
xlim([0.5 num_sub+0.5]);
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
ax.XTick = linspace(1,num_sub,num_sub);
hold on;

subplot(2,3,3)
for i=1:(num_sub)
    if sum(i==authors)==1
        plot(1:frames,squeeze(temporal_kernel(i,1:frames)),'--k');
    else
        plot(1:frames,squeeze(temporal_kernel(i,1:frames)),'k');
    end
    hold on;
end
LH(1) = plot(1:frames,mean(squeeze(temporal_kernel(:,1:frames)),1),'-ok','LineWidth',2);
L{1} = 'Foveated';
xlabel('Fixation Number','Fontsize',20);
xticks([1 2 3 4]);
ylabel({'Weights'},'Fontsize',20);
hold on;
yline(0.0,'k','LineWidth',1.5);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
hold on;
for i=1:(num_sub)
    if sum(i==authors)==1
        plot((1:frames),[squeeze(temporal_kernel(i,(frames+1):(end - pr ))) mean(squeeze(temporal_kernel(i,(end - pr+1 ):end)))],'-dm');
        hold on;
        scatter([frames frames]',squeeze(temporal_kernel(i,(end - pr+1 ):end)),'dm');
    else
        plot((1:frames),[squeeze(temporal_kernel(i,(frames+1):(end - pr ))) mean(squeeze(temporal_kernel(i,(end - pr+1 ):end)))],'-om');
        hold on;
        scatter([frames frames]',squeeze(temporal_kernel(i,(end - pr+1 ):end)),'om');
    end
end
temp1 = mean(squeeze(temporal_kernel(:,(frames+1):(end - pr ))),1);
temp2 = mean(mean(squeeze(temporal_kernel(:,(end - pr+1 ):end)),1));
temp = [temp1 temp2];
LH(2) = plot((1:frames),temp,'-om','LineWidth',2);
L{2} = 'Non-foveated';
scatter([frames frames]',mean(squeeze(temporal_kernel(:,(end - pr+1 ):end)),1),100,'om','filled');
xlim([1 4]);
ylim([-2 5]);
legend(LH,L, 'Fontsize',20, 'Box','off');

sgtitle('All subject all analysis w.r.t physical time','fontsize',30);

%%
figure();
subplot(2,1,1)
k1 = 1;
k2 = 2;
k = [];
for sub=1:num_sub
    if sub==num_sub
        LH6(1) = bar(k1, squeeze(sub_accuracy_first_half(sub)),'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',0.5);
        L6{1} = 'First half';
    else
        if sum(sub==chosen_sub)==1
            bar(k1, squeeze(sub_accuracy_first_half(sub)),'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',3);
        elseif sum(sub==authors)==1
            bar(k1, squeeze(sub_accuracy_first_half(sub)),'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor',[0.6350 0.0780 0.1840],'LineWidth',3);
        else
            bar(k1, squeeze(sub_accuracy_first_half(sub)),'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',0.5);
        end
    end
    hold on;
    if sub==num_sub
        LH6(2) = bar(k2, sub_accuracy_second_half(sub),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.5);
        L6{2} = 'Second half';
    else
        if sum(sub==chosen_sub)==1
            bar(k2, sub_accuracy_second_half(sub),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',3);
        elseif sum(sub==authors)==1
            bar(k2, sub_accuracy_second_half(sub),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor',[0.6350 0.0780 0.1840],'LineWidth',3);
        else
            bar(k2, sub_accuracy_second_half(sub),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.5);
        end
    end
    hold on;
    err1 = sqrt((sub_accuracy_first_half(sub)*(1-sub_accuracy_first_half(sub)))/floor(trials_per_num_peri(sub)/2));
    err2 = sqrt((sub_accuracy_second_half(sub)*(1-sub_accuracy_second_half(sub)))/(trials_per_num_peri(sub) -floor(trials_per_num_peri(sub)/2)));
    errorbar([k1 k2],[squeeze(sub_accuracy_first_half(sub)) sub_accuracy_second_half(sub)],[err1 err2],'ok','LineWidth',2,'linestyle','none');
    xlabel('Subject number','Fontsize',20);
    ylabel({'Percent Correct'},'Fontsize',20);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.4 1]);
    yline(0.5,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    k = [k; k1];
    k1 = k1 + 2.5;
    k2 = k2 + 2.5;
end
xticklabels(linspace(1, num_sub, num_sub));
xticks(k+0.5);
legend(LH6,L6, 'Fontsize',20, 'Box','off');
hold on;

subplot(2,1,2)
k1 = 1;
k2 = 2;
k = [];
for sub=1:num_sub
    if sub==num_sub
        LH6(1) = bar(k1, squeeze(prob_chose_in_favor_opp_first_half(sub)),'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',0.5);
        L6{1} = 'First half';
    else
        if sum(sub==chosen_sub)==1
            bar(k1, squeeze(prob_chose_in_favor_opp_first_half(sub)),'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',3);
        elseif sum(sub==authors)==1
            bar(k1, squeeze(prob_chose_in_favor_opp_first_half(sub)),'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor',[0.6350 0.0780 0.1840],'LineWidth',3);
        else
            bar(k1, squeeze(prob_chose_in_favor_opp_first_half(sub)),'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',0.5);
        end
    end
    hold on;
    if sub==num_sub
        LH6(2) = bar(k2, prob_chose_in_favor_opp_second_half(sub),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.5);
        L6{2} = 'Second half';
    else
        if sum(sub==chosen_sub)==1
            bar(k2, prob_chose_in_favor_opp_second_half(sub),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',3);
        elseif sum(sub==authors)==1
            bar(k2, prob_chose_in_favor_opp_second_half(sub),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor',[0.6350 0.0780 0.1840],'LineWidth',3);
        else
            bar(k2, prob_chose_in_favor_opp_second_half(sub),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.5);
        end
    end
    hold on;
    errorbar([k1 k2],[squeeze(prob_chose_in_favor_opp_first_half(sub)) prob_chose_in_favor_opp_second_half(sub)],[squeeze(lo_err_opp_first_half(sub)) squeeze(lo_err_opp_second_half(sub))],[squeeze(hi_err_opp_first_half(sub)) squeeze(hi_err_opp_second_half(sub))],'ok','LineWidth',2,'linestyle','none');
    xlabel('Subject number','Fontsize',20);
    ylabel({'Probability of saccading to ';'confirming evidence'},'Fontsize',20);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.4 1]);
    yline(0.5,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    k = [k; k1];
    k1 = k1 + 2.5;
    k2 = k2 + 2.5;
end
xticklabels(linspace(1, num_sub, num_sub));
xticks(k+0.5);
legend(LH6,L6, 'Fontsize',20, 'Box','off');

sgtitle('All subject all analysis w.r.t first half and second half trials','fontsize',30);

%%
figure();
subplot(2,3,1)
hold on;
bar(1:num_sub,sub_accuracy./trials_per_num_peri, 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',0.75);
hold on;
bar(authors,sub_accuracy(authors)./trials_per_num_peri(authors), 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor',[0.6350 0.0780 0.1840],'LineWidth',3.0);
hold on;
yline(0.5,'k','LineWidth',1.5);
xlabel('Subject Number','Fontsize',20);
ylabel('Percent Correct','Fontsize',20);
ax = gca;
ax.LineWidth=2;
ylim([0.4 1]);
xlim([0.5 num_sub+0.5]);
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
ax.XTick = linspace(1,num_sub,num_sub);
hold on;
tmp = sub_accuracy./trials_per_num_peri;
errorbar(1:num_sub,tmp,sqrt((tmp.*(1-tmp))./trials_per_num_peri),'ok','LineWidth', 1.5);
hold on;

subplot(2,3,2)
hold on;
bar(1:num_sub,squeeze(prob_chose_in_favor_opp{1}), 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',0.75);
bar(authors,squeeze(prob_chose_in_favor_opp{1}([2 3])), 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor',[0.6350 0.0780 0.1840],'LineWidth',3);
hold on;
errorbar(1:num_sub,squeeze(prob_chose_in_favor_opp{1}),squeeze(lo_err_opp{1}),squeeze(hi_err_opp{1}),'ok','LineWidth', 1.5)
hold on;
text([1:num_sub]-0.3,squeeze(prob_chose_in_favor_opp{1})+0.065,stars,'FontSize',15,'FontWeight','bold');
xlabel('Subject Number','Fontsize',20)
ylabel({'Probability of saccading to ';'confirming evidence'},'Fontsize',20);
yline(0.5,'k','LineWidth',1.5);
ax = gca;
ax.LineWidth=2;
ylim([0.4 1]);
xlim([0.5 num_sub+0.5]);
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
ax.XTick = linspace(1,num_sub,num_sub);
hold on;

subplot(2,3,3)
for i=1:(num_sub)
    if sum(i==authors)==1
        plot(1:frames,squeeze(temporal_kernel(i,1:frames)),'--k');
    else
        plot(1:frames,squeeze(temporal_kernel(i,1:frames)),'k');
    end
    hold on;
end
LH(1) = plot(1:frames,mean(squeeze(temporal_kernel(:,1:frames)),1),'-ok','LineWidth',2);
L{1} = 'Foveated';
xlabel('Fixation Number','Fontsize',20);
xticks([1 2 3 4]);
ylabel({'Weights'},'Fontsize',20);
hold on;
yline(0.0,'k','LineWidth',1.5);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
hold on;
for i=1:(num_sub)
    if sum(i==authors)==1
        scatter((2:frames),squeeze(temporal_kernel(i,(frames+1):(end - pr ))),'dm');
        hold on;
        scatter([frames+1 frames+1]',squeeze(temporal_kernel(i,(end - pr+1 ):end)),'dm');
    else
        scatter((2:frames),squeeze(temporal_kernel(i,(frames+1):(end - pr ))),'om');
        hold on;
        scatter([frames+1 frames+1]',squeeze(temporal_kernel(i,(end - pr+1 ):end)),'om');
    end
end
LH(2) = scatter((2:frames),mean(squeeze(temporal_kernel(:,(frames+1):(end - pr ))),1),100,'m','filled');
L{2} = 'Non-foveated';
scatter([frames+1 frames+1]',mean(squeeze(temporal_kernel(:,(end - pr+1 ):end)),1),100,'m','filled');
xlim([1 5]);
ylim([-1 5]);
legend(LH,L, 'Fontsize',20, 'Box','off');

sgtitle('All subject all analysis','fontsize',30);
%%
figure();
k1 = 1;
k2 = 2;
k = [];
for sub=1:num_sub
    if sub==num_sub
        LH6(1) = bar(k1, squeeze(prob_chose_in_favor_opp{1}(sub,:)),'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',0.5);
        L6{1} = 'Bias w.r.t. accumulated belief';
    else
        if sum(sub==chosen_sub)==1
            bar(k1, squeeze(prob_chose_in_favor_opp{1}(sub,:)),'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',3);
        elseif sum(sub==authors)==1
            bar(k1, squeeze(prob_chose_in_favor_opp{1}(sub,:)),'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor',[0.6350 0.0780 0.1840],'LineWidth',3);
        else
            bar(k1, squeeze(prob_chose_in_favor_opp{1}(sub,:)),'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',0.5);
        end
    end
    hold on;
    if sub==num_sub
        LH6(2) = bar(k2, opp_raw_prob{1}(sub,:),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.5);
        L6{2} = 'Bias w.r.t foveated stimuli only';
    else
        if sum(sub==chosen_sub)==1
            bar(k2, opp_raw_prob{1}(sub,:),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',3);
        elseif sum(sub==authors)==1
            bar(k2, opp_raw_prob{1}(sub,:),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor',[0.6350 0.0780 0.1840],'LineWidth',3);
        else
            bar(k2, opp_raw_prob{1}(sub,:),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.5);
        end
    end
    hold on;
    errorbar([k1 k2],[squeeze(prob_chose_in_favor_opp{1}(sub,:)) opp_raw_prob{1}(sub,:)],[squeeze(lo_err_opp{1}(sub,:)) opp_raw_lo_err{1}(sub,:)],[squeeze(hi_err_opp{1}(sub,:)) opp_raw_hi_err{1}(sub,:)],'ok','LineWidth',2,'linestyle','none');
    xlabel('Subject number','Fontsize',20);
    ylabel({'Probability of saccading to ';'confirming evidence'},'Fontsize',20);
    hold on;
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ylim([0.4 1]);
    yline(0.5,'-k','linewidth',2);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    k = [k; k1];
    k1 = k1 + 2.5;
    k2 = k2 + 2.5;
end
xticklabels(linspace(1, num_sub, num_sub));
xticks(k+0.5);
legend(LH6,L6, 'Fontsize',20, 'Box','off');
title('Pattern Matching: Bias w.r.t belief vs foveated stimuli','fontsize',30);

%%
figure();
for i=1:(num_sub)
    if sum(i==chosen_sub)==1
        bar(i,squeeze(-pattern_one_previous{i}(4,1) + belief_one_previous{i}(4,1)),'FaceColor',[0.75, 0.75, 0],'EdgeColor','k','LineWidth',3);
        hold on;
    elseif sum(i==authors)==1
        bar(i,squeeze(-pattern_one_previous{i}(4,1) + belief_one_previous{i}(4,1)),'FaceColor',[0.75, 0.75, 0],'EdgeColor',[0.6350 0.0780 0.1840],'LineWidth',3);
        hold on;
    else
        bar(i,squeeze(-pattern_one_previous{i}(4,1) + belief_one_previous{i}(4,1)),'FaceColor',[0.75, 0.75, 0],'EdgeColor','k','LineWidth',0.75);
        hold on;
    end
    errorbar(i,squeeze(-pattern_one_previous{i}(4,1) + belief_one_previous{i}(4,1)),sqrt(squeeze(pattern_one_previous{i}(4,2)^2 + belief_one_previous{i}(4,2))^2),'-ko','LineWidth',2,'linestyle','none');
    hold on;
    yline(0.0,'k','LineWidth',2);
    ylabel('Probability (belief based - current stimulus)','fontsize',30);
    xlabel('Subject number','Fontsize',30);
    xlim([0.25 num_sub+0.75]);
    ylim([-0.2 0.15]);
    xticks(linspace(1,num_sub,num_sub));
    ax = gca;
    set(ax, 'box','off');
    ax.LineWidth=2;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
end
sgtitle('Pattern Matching: Bias w.r.t current belief vs current foveated stimulus','fontsize',30);

% figure();
% for i=1:(num_sub)
%     if sum(i==chosen_sub)==1
%         bar(i,squeeze(-pattern_one_previous_case_belief_evi_opp{i}(4,1) + belief_one_previous_case_belief_evi_opp{i}(4,1)),'FaceColor',[0.75, 0.75, 0],'EdgeColor','k','LineWidth',3);
%         hold on;
%     elseif sum(i==authors)==1
%         bar(i,squeeze(-pattern_one_previous_case_belief_evi_opp{i}(4,1) + belief_one_previous_case_belief_evi_opp{i}(4,1)),'FaceColor',[0.75, 0.75, 0],'EdgeColor',[0.6350 0.0780 0.1840],'LineWidth',3);
%         hold on;
%     else
%         bar(i,squeeze(-pattern_one_previous_case_belief_evi_opp{i}(4,1) + belief_one_previous_case_belief_evi_opp{i}(4,1)),'FaceColor',[0.75, 0.75, 0],'EdgeColor','k','LineWidth',0.75);
%         hold on;
%     end
%     errorbar(i,squeeze(-pattern_one_previous_case_belief_evi_opp{i}(4,1) + belief_one_previous_case_belief_evi_opp{i}(4,1)),sqrt(squeeze(pattern_one_previous_case_belief_evi_opp{i}(4,2)^2 + belief_one_previous_case_belief_evi_opp{i}(4,2))^2),'-ko','LineWidth',2,'linestyle','none');
%     hold on;
%     yline(0.0,'k','LineWidth',2);
%     ylabel('Probability (belief based - current stimulus)','fontsize',30);
%     xlabel('Subject number','Fontsize',30);
%     xlim([0.25 num_sub+0.75]);
%     ylim([-0.15 0.15]);
%     xticks(linspace(1,num_sub,num_sub));
%     ax = gca;
%     set(ax, 'box','off');
%     ax.LineWidth=2;
%     ax.XAxis.FontSize = 20;
%     ax.YAxis.FontSize = 20;
% end
% sgtitle({'Pattern Matching when belief and current stimulus have opposite sign:';'Bias w.r.t current belief vs current foveated stimulus'},'fontsize',30);

figure();
for i=1:(num_sub)
    if (i==chosen_sub)==1
        bar(i,squeeze(-opp_raw_prob{1}(i,:)) + squeeze(prob_chose_in_favor_opp{1}(i,:)),'FaceColor',[0.75, 0.75, 0],'EdgeColor','k','LineWidth',3);
        hold on;
    elseif sum(i==authors)==1
        bar(i,squeeze(-opp_raw_prob{1}(i,:)) + squeeze(prob_chose_in_favor_opp{1}(i,:)),'FaceColor',[0.75, 0.75, 0],'EdgeColor',[0.6350 0.0780 0.1840],'LineWidth',3);
        hold on;
    else
        bar(i,squeeze(-opp_raw_prob{1}(i,:)) + squeeze(prob_chose_in_favor_opp{1}(i,:)),'FaceColor',[0.75, 0.75, 0],'EdgeColor','k','LineWidth',0.75);
        hold on;
    end
    hold on;
    yline(0.0,'k','LineWidth',2);
    ylabel('Probability (belief based - current stimulus)','fontsize',30);
    xlabel('Subject number','Fontsize',30);
    xlim([0.25 num_sub+0.75]);
    ylim([-0.2 0.2]);
    xticks(linspace(1,num_sub,num_sub));
    ax = gca;
    set(ax, 'box','off');
    ax.LineWidth=2;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
end
errorbar(1:num_sub,squeeze(-opp_raw_prob{1}) + squeeze(prob_chose_in_favor_opp{1}),sqrt(squeeze((lo_err_opp{1}).^2 + (opp_raw_lo_err{1}).^2)),'-ko','LineWidth',2,'linestyle','none');
% errorbar(1:num_sub,squeeze(-opp_raw_prob{1}) + squeeze(prob_chose_in_favor_opp{1}),sqrt(squeeze((lo_err_opp{1}+hi_err_opp{1}).^2 + (opp_raw_lo_err{1}+opp_raw_hi_err{1}).^2)),'-ko','LineWidth',2,'linestyle','none');
sgtitle('Pattern Matching: Bias w.r.t current belief vs all foveated stimuli so far','fontsize',30);

%%
figure();
subplot(1,2,1)
for i=1:(num_sub)
    if sum(i==authors)==3
        plot([1 2],squeeze(prob_chose_in_favor_opp{2}(i,:)),'-ro','LineWidth',1.0);
    else
        plot([1 2],squeeze(prob_chose_in_favor_opp{2}(i,:)),'-bo','LineWidth',1.0);
    end
    hold on;
end
hold on;
scatter(ones(1,num_sub),squeeze(prob_chose_in_favor_opp{2}(:,1)),100,'b','filled');
hold on;
scatter(ones(1,2),squeeze(prob_chose_in_favor_opp{2}(authors,1)),100,'r','filled');
hold on;
LH(1) = scatter(ones(1,num_sub)*2,squeeze(prob_chose_in_favor_opp{2}(:,2)),100,'b','filled');
L{1} = 'Naive Subject';
hold on;
LH(2) = scatter(ones(1,2)*2,squeeze(prob_chose_in_favor_opp{2}(authors,2)),100,'r','filled');
L{2} = 'Author';
hold on;
plot([1 2],squeeze(mean(prob_chose_in_favor_opp{2},1)),'-ko','LineWidth',5);
yline(0.5,'k','LineWidth',2);
xlabel('Accumulated evidence','Fontsize',20);
ylabel('Prob chose in favor','Fontsize',20);
xlim([0.75 2.25]);
ylim([0.3 1.0]);
set(gca, 'XTick', [1 2]);
set(gca, 'XTickLabel',{'Less evidence' 'More evidence'});
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
hold on;
legend(LH,L, 'Fontsize',20, 'Box','off');

subplot(1,2,2)
errorbar(squeeze(prob_chose_in_favor_opp{2}(:,1)),squeeze(prob_chose_in_favor_opp{2}(:,2)),squeeze(lo_err_opp{2}(:,1)),squeeze(hi_err_opp{2}(:,1)),'horizontal','LineWidth',2, 'LineStyle', 'none','Color','k');
hold on;
errorbar(squeeze(prob_chose_in_favor_opp{2}(:,1)),squeeze(prob_chose_in_favor_opp{2}(:,2)),squeeze(lo_err_opp{2}(:,2)),squeeze(hi_err_opp{2}(:,2)),'vertical', 'LineWidth',2,'LineStyle', 'none','Color','k');
hold on;
scatter(squeeze(prob_chose_in_favor_opp{2}(:,1)),squeeze(prob_chose_in_favor_opp{2}(:,2)),100,'b','filled');
hold on;
scatter(squeeze(prob_chose_in_favor_opp{2}([2 3],1)),squeeze(prob_chose_in_favor_opp{2}([2 3],2)),100,'r','filled');
hold on;
plot(linspace(0.3,1.0,100),linspace(0.3,1.0,100),'k','LineWidth',2);
xlabel({'Probability of biased saccade';' for low evidence accumulated'},'Fontsize',20);
ylabel({'Probability of biased saccade';' for high evidence accumulated'},'Fontsize',20);
xlim([0.3 1.0]);
ylim([0.3 1.0]);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;

%%

figure();
subplot(1,2,1)
cnt = 1;
hold on;
LH(1) = scatter(mean(temporal_kernel(:,8:9),2),squeeze(prob_chose_in_favor_opp{1}),100,'b','filled');
hold on;
errorbar(mean(temporal_kernel(:,8:9),2),squeeze(prob_chose_in_favor_opp{1}),squeeze(lo_err_opp{1}),squeeze(hi_err_opp{1}),'vertical','LineWidth',2, 'LineStyle', 'none','Color','k');
hold on;
errorbar(mean(temporal_kernel(:,8:9),2),squeeze(prob_chose_in_favor_opp{1}),squeeze(std(temporal_kernel(:,8:9)'))/sqrt(2),'horizontal','LineWidth',2, 'LineStyle', 'none','Color','k');
hold on;
LH(2) = scatter(mean(temporal_kernel(chosen_sub,8:9),2),squeeze(prob_chose_in_favor_opp{1}(chosen_sub)),100,'r','filled');
hold on;
errorbar(mean(temporal_kernel(chosen_sub,8:9),2),squeeze(prob_chose_in_favor_opp{1}(chosen_sub)),squeeze(lo_err_opp{1}(chosen_sub)),squeeze(hi_err_opp{1}(chosen_sub)),'vertical','LineWidth',2, 'LineStyle', 'none','Color','k');
hold on;
errorbar(mean(temporal_kernel(chosen_sub,8:9),2),squeeze(prob_chose_in_favor_opp{1}(chosen_sub)),squeeze(std(temporal_kernel(chosen_sub,8:9)'))/sqrt(2),'horizontal','LineWidth',2, 'LineStyle', 'none','Color','k');
hold on;
L{1} = 'Biased Subject';
hold on;
L{2} = 'Not biased Subject';
hold on;
yline(0.5,'k','LineWidth',2);
xlabel('Mean of last two peripheral weights','Fontsize',20);
ylabel('Prob chose in favor','Fontsize',20);
[r,p] = corrcoef(mean(temporal_kernel(:,8:9),2),squeeze(prob_chose_in_favor_opp{1}));
text(0.75,0.75,['corr=' num2str(r(1,2))],'fontsize',15);
text(0.75,0.7,['p-val=' num2str(p(1,2))],'fontsize',15);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
hold on;
legend(LH,L, 'Fontsize',20, 'Box','off');


subplot(1,2,2)
cnt = 1;
hold on;
LH(1) = scatter(mean(temporal_kernel(:,5:end),2),squeeze(prob_chose_in_favor_opp{1}),100,'b','filled');
hold on;
errorbar(mean(temporal_kernel(:,5:end),2),squeeze(prob_chose_in_favor_opp{1}),squeeze(lo_err_opp{1}),squeeze(hi_err_opp{1}),'vertical','LineWidth',2, 'LineStyle', 'none','Color','k');
hold on;
errorbar(mean(temporal_kernel(:,5:end),2),squeeze(prob_chose_in_favor_opp{1}),squeeze(std(temporal_kernel(:,5:end)'))/sqrt(5),'horizontal','LineWidth',2, 'LineStyle', 'none','Color','k');
hold on;
LH(2) = scatter(mean(temporal_kernel(chosen_sub,5:end),2),squeeze(prob_chose_in_favor_opp{1}(chosen_sub)),100,'r','filled');
hold on;
errorbar(mean(temporal_kernel(chosen_sub,5:end),2),squeeze(prob_chose_in_favor_opp{1}(chosen_sub)),squeeze(lo_err_opp{1}(chosen_sub)),squeeze(hi_err_opp{1}(chosen_sub)),'vertical','LineWidth',2, 'LineStyle', 'none','Color','k');
hold on;
errorbar(mean(temporal_kernel(chosen_sub,5:end),2),squeeze(prob_chose_in_favor_opp{1}(chosen_sub)),squeeze(std(temporal_kernel(chosen_sub,5:end)'))/sqrt(5),'horizontal','LineWidth',2, 'LineStyle', 'none','Color','k');
hold on;
L{1} = 'Biased Subject';
hold on;
L{2} = 'Not biased Subject';
hold on;
yline(0.5,'k','LineWidth',2);
xlabel('Mean of all peripheral weights','Fontsize',20);
ylabel('Prob chose in favor','Fontsize',20);
[r,p] = corrcoef(mean(temporal_kernel(:,5:end),2),squeeze(prob_chose_in_favor_opp{1}));
text(0.5,0.75,['corr=' num2str(r(1,2))],'fontsize',15);
text(0.5,0.7,['p-val=' num2str(p(1,2))],'fontsize',15);
ax = gca;
ax.LineWidth=2;
set(ax, 'box','off');
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
hold on;
legend(LH,L, 'Fontsize',20, 'Box','off');
%%
chosen_sub = [];
for i=1:num_sub
    if stars{i}==' '
        chosen_sub(end+1) = i;
    end
end
if ~isempty(chosen_sub)
    figure();
    subplot(2,3,1)
    hold on;
    bar(1:length(chosen_sub),sub_accuracy(chosen_sub)./trials_per_num_peri(chosen_sub), 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',0.75);
    hold on;
    yline(0.5,'k','LineWidth',1.5);
    xlabel('Subject Number','Fontsize',20);
    ylabel('Percent Correct','Fontsize',20);
    ax = gca;
    ax.LineWidth=2;
    ylim([0.4 1]);
    xlim([0.5 length(chosen_sub)+0.5]);
    ax.YAxis.FontSize = 20;
    ax.XAxis.FontSize = 20;
    xticks(linspace(1, length(chosen_sub),length(chosen_sub)));
    xticklabels(chosen_sub);
    hold on;
    tmp = sub_accuracy(chosen_sub)./(trials_per_num_peri(chosen_sub));
    errorbar(1:length(chosen_sub),tmp,tmp.*(1-tmp)./2,'ok','LineWidth', 1.5);
    hold on;
    
    subplot(2,3,2)
    hold on;
    bar(1:length(chosen_sub),squeeze(squeeze(prob_chose_in_favor_opp{1}(chosen_sub,:))), 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',0.75);
    hold on;
    errorbar(1:length(chosen_sub),squeeze(squeeze(prob_chose_in_favor_opp{1}(chosen_sub,:))),squeeze(lo_err_opp{1}(chosen_sub,:)),squeeze(hi_err_opp{1}(chosen_sub,:)),'ok','LineWidth', 1.5);
    hold on;
    xlabel('Subject Number','Fontsize',20);
    ylabel({'Probability of saccading to ';'evidence confirming current belief'},'Fontsize',20);
    yline(0.5,'k','LineWidth',1.5);
    ax = gca;
    ax.LineWidth=2;
    ylim([0.4 1]);
    xticks(linspace(1, length(chosen_sub),length(chosen_sub)));
    xticklabels(chosen_sub);
    xlim([0.5 length(chosen_sub)+0.5]);
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    hold on;
    
    subplot(2,3,3)
    for i=1:length(chosen_sub)
        plot(1:frames,squeeze(temporal_kernel(chosen_sub(i),1:frames)),'k','linewidth',2);
        hold on;
    end
    % ylim([-0.05 0.05]);
    xlabel('Fixation Number','Fontsize',20);
    xticks([1 2 3 4]);
    ylabel({'Weights given to stimuli'; 'in final decision'},'Fontsize',20);
    hold on;
    yline(0.0,'k','LineWidth',1.5);
    ax = gca;
    ax.LineWidth=2;
    set(ax, 'box','off');
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    hold on;
    for i=1:length(chosen_sub)
        scatter((2:frames),squeeze(temporal_kernel(chosen_sub(i),(frames+1):(end - pr ))),'om','filled');
        hold on;
        scatter([frames+1 frames+1]',squeeze(temporal_kernel(chosen_sub(i),(end - pr+1 ):end)),'om','filled');
    end
    xlim([1 5]);
    
    sgtitle('Analysis for subjects with no effect','fontsize',30);
end
%%
figure();
for i=1:num_sub
    subplot(2,sb_plt,i)
    errorbar(mn_sig_ideal_all_norm(i,:),mn_perf_ideal_all(i,:),err_perf_ideal_all(i,:),'-og','linewidth',2);
    hold on;
    errorbar(mn_sig_norm(i,:),mn_perf(i,:),err_perf(i,:),'-ob','linewidth',2);
    hold on;
    errorbar(mn_sig_ideal_foveal_only_norm(i,:),mn_perf_ideal_foveal_only(i,:),err_perf_ideal_foveal_only(i,:),'-or','linewidth',2);
    hold on;
    title(['Subject ' num2str(i)],'fontsize',20);
    if i==md_x || i==(sb_plt+md_x)
        xlabel('Normalized Evidence Accumulated','fontsize',20);
    end
    if i==1 || i==(sb_plt+1)
        ylabel('Probability of correct answer','fontsize',20);
    end
    ylim([0.3 1.001]);
    xlim([0.0 1.001]);
    yline(0.5,'k','Linewidth',2);
    ax = gca;
    set(ax, 'box','off');
    ax.LineWidth=2;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
end

sgtitle(['Green: Ideal all frames Red: Ideal fovea only Blue: Subject'], 'fontsize',30);

%%
figure();
for i=1:num_sub
    subplot(2,sb_plt,i)
    errorbar(mn_sig_ideal_all_actual_sig(i,:),mn_perf_ideal_all_actual_sig(i,:),err_perf_ideal_all_actual_sig(i,:),'-og','linewidth',2);
    hold on;
    errorbar(mn_sig_ideal_all_actual_sig(i,:),mn_perf_actual_sig(i,:),err_perf_actual_sig(i,:),'-ob','linewidth',2);
    hold on;
    errorbar(mn_sig_ideal_all_actual_sig(i,:),mn_perf_ideal_foveal_only_actual_sig(i,:),err_perf_ideal_foveal_only_actual_sig(i,:),'-or','linewidth',2);
    hold on;
    title(['Subject ' num2str(i)],'fontsize',20);
    if i==md_x || i==(sb_plt+md_x)
        xlabel('Raw evidence shown','fontsize',20);
    end
    if i==1 || i==(sb_plt+1)
        ylabel('Probability of correct answer','fontsize',20);
    end
    ylim([0.45 1.001]);
    yline(0.5,'k','Linewidth',2);
    xticks(mn_sig_ideal_all_actual_sig(i,:));
    set(gca,'xticklabel',num2str(get(gca,'xtick')','%.1f'))
    ax = gca;
    set(ax, 'box','off');
    ax.LineWidth=2;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
end
sgtitle(['Green: Ideal all frames Red: Ideal fovea only Blue: Subject'], 'fontsize',30);

%%
% bin_num = 3;
% figure();
% llo_mn = [];
% choice_mn = [];
% err_ch = [];
% err_ch_mn = [];
% llo_mean = [];
% choice_mean = [];
% for sub=1:num_sub
%     subplot(2,7,sub)
%     %     llo = llo( ~any( isnan( llo ) | isinf( llo ), 2 ),: );
%     [sorted_llo,order_llo] = sort(log_bernoulli{sub});
%     choice_used = choice{sub}(order_llo);
%     bin_edges = linspace(min(sorted_llo),max(sorted_llo),bin_num+1);
%     for bn=1:length(bin_edges)-1
%         llo_mean(bn) = mean(sorted_llo(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)));
%         choice_mean(bn) = mean(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)));
%         err_ch(bn) = sqrt(var(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1)))/length(choice_used(sorted_llo>=bin_edges(bn) & sorted_llo<=bin_edges(bn+1))));
%     end
%     [llo_mn, ord] = sort(llo_mean);
%     choice_mn = choice_mean(ord);
%     err_ch_mn = err_ch(ord);
%     errorbar(llo_mn,choice_mn,err_ch_mn,'-or','Linewidth',2);
%     %     shadedErrorBar(llo_mn,choice_mn,err_ch_mn);
%     if sub==1 || sub==6
%         ylabel('Probability chose vertical','Fontsize',20);
%     end
%     if sub==3 || sub==8
%         xlabel('Log likelihood odds','Fontsize',20);
%     end
%     hold on;
%     ax = gca;
%     ax.LineWidth=2;
%     set(ax, 'box','off');
%     ylim([0.0 1]);
%     xlim([-1*max(abs(llo_mn)) max(abs(llo_mn))]);
%     yline(0.5,'-k','linewidth',2);
%     hold on;
%     xline(0.0,'-k','linewidth',2);
%     ax.XAxis.FontSize = 20;
%     ax.YAxis.FontSize = 20;
%     title(['Subject ' num2str(sub)], 'fontsize',20)
% end
% sgtitle('Check how good weights predict behavior','fontsize',30);

%%
% figure();
% for i=1:(num_sub)
%     s=subplot(2,7,i);
%     for cs = 1:4
%         bar(cs,squeeze(pattern_one_previous{i}(cs,1)),'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.75);
%         hold on;
%         errorbar(cs,squeeze(pattern_one_previous{i}(cs,1)),squeeze(pattern_one_previous{i}(cs,2)),'-ko','LineWidth',2,'linestyle','none');
%         hold on;
%         yline(0.5,'k','LineWidth',2);
%         if i==1 || i==6
%             ylabel('Probability','fontsize',30)
%         end
%     end
%     xlim([0.25 4.75]);
%     ylim([0.3 1.0]);
%     xticks([1 2 3 4]);
%     xticklabels({'1sc', '2sc', '3sc', 'all'});
%     ax = gca;
%     set(ax, 'box','off');
%     ax.LineWidth=2;
%     ax.XAxis.FontSize = 10;
%     ax.YAxis.FontSize = 20;
%     title(['Subject ' num2str(i)],'fontsize',20);
% end
% sgtitle('One previous pattern matching','fontsize',30);
% %
% % %%
% figure();
% for i=1:(num_sub)
%     s=subplot(2,7,i);
%     for cs = 1:4
%         bar(cs,squeeze(pattern_one_previous{i}(cs,1) - belief_one_previous{i}(cs,1)),'FaceColor',[0.75, 0.75, 0],'EdgeColor','k','LineWidth',0.75);
%         hold on;
%         errorbar(cs,squeeze(pattern_one_previous{i}(cs,1) - belief_one_previous{i}(cs,1)),sqrt(squeeze(pattern_one_previous{i}(cs,2)^2 + belief_one_previous{i}(cs,2))^2),'-ko','LineWidth',2,'linestyle','none');
%         hold on;
%         yline(0.0,'k','LineWidth',2);
%         if i==1 || i==6
%             ylabel('Probability','fontsize',30)
%         end
%     end
%     xlim([0.25 4.75]);
%     ylim([-1.0 1.0]);
%     xticks([1 2 3 4]);
%     xticklabels({'1sc', '2sc', '3sc', 'all'});
%     ax = gca;
%     set(ax, 'box','off');
%     ax.LineWidth=2;
%     ax.XAxis.FontSize = 10;
%     ax.YAxis.FontSize = 20;
%     title(['Subject ' num2str(i)],'fontsize',20);
% end
% sgtitle('One previous difference between belief and pattern matching','fontsize',30);
% 
% %%
% figure();
% for i=1:(num_sub)
%     s=subplot(2,7,i);
%     for cs = 1:3
%         bar(cs,squeeze(pattern_three_previous{i}(cs,1)),'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.75);
%         hold on;
%         errorbar(cs,squeeze(pattern_three_previous{i}(cs,1)),squeeze(pattern_three_previous{i}(cs,2)),'-ko','LineWidth',2,'linestyle','none');
%         hold on;
%         yline(0.5,'k','LineWidth',2);
%         if i==1 || i==6
%             ylabel('Probability','fontsize',30)
%         end
%     end
%     xlim([0.25 3.75]);
%     ylim([0.3 1.0]);
%     xticks([1 2 3]);
%     xticklabels({'3sc 3sm', '3sc 2sm', 'all'});
%     ax = gca;
%     set(ax, 'box','off');
%     ax.LineWidth=2;
%     ax.XAxis.FontSize = 10;
%     ax.YAxis.FontSize = 20;
%     title(['Subject ' num2str(i)],'fontsize',20);
% end
% sgtitle('Three previous pattern matching','fontsize',30);
%
% %%
% figure();
% for i=1:(num_sub)
%     s=subplot(2,7,i);
%     for cs = 1:3
%         bar(cs,squeeze(pattern_three_previous{i}(cs,1) - belief_three_previous{i}(cs,1)),'FaceColor',[0.75, 0.75, 0],'EdgeColor','k','LineWidth',0.75);
%         hold on;
%         errorbar(cs,squeeze(pattern_three_previous{i}(cs,1) - belief_three_previous{i}(cs,1)),sqrt(squeeze(pattern_three_previous{i}(cs,2)^2 + belief_three_previous{i}(cs,2))^2),'-ko','LineWidth',2,'linestyle','none');
%         hold on;
%         yline(0.0,'k','LineWidth',2);
%         if i==1 || i==6
%             ylabel('Probability','fontsize',30)
%         end
%     end
%     xlim([0.25 3.75]);
%     ylim([-1.0 1.0]);
%     xticks([1 2 3]);
%     xticklabels({'3sc 3sm', '3sc 2sm', 'all'});
%     ax = gca;
%     set(ax, 'box','off');
%     ax.LineWidth=2;
%     ax.XAxis.FontSize = 10;
%     ax.YAxis.FontSize = 20;
%     title(['Subject ' num2str(i)],'fontsize',20);
% end
% sgtitle('Three previous difference between belief and pattern matching','fontsize',30);
%
% %%
% figure();
% for i=1:(num_sub)
%     s=subplot(2,7,i);
%     cnt = 0;
%     for cs = 1:3
%         for cs_cs = 1:3
%             cnt = cnt + 1;
%             bar(cnt,squeeze(pattern_two_previous{i}(cs,cs_cs,1)),'FaceColor',[0, 0.4470, 0.7410],'EdgeColor','k','LineWidth',0.75);
%             hold on;
%             errorbar(cnt,squeeze(pattern_two_previous{i}(cs,cs_cs,1)),squeeze(pattern_two_previous{i}(cs,cs_cs,2)),'-ko','LineWidth',2,'linestyle','none');
%             hold on;
%             yline(0.5,'k','LineWidth',2);
%             if i==1 || i==6
%                 ylabel('Probability','fontsize',30)
%             end
%         end
%     end
%     xlim([0.25 9.75]);
%     ylim([0.0 1.0]);
%     xticks([1 2 3 4 5 6 7 8 9]);
%     row1 = {'2sc', '2sc', '2sc',...
%         '3sc', '3sc', '3sc'...
%         'all', 'all', 'all'};
%     row2 = {'sm', 'df', 'all',...
%         'sm', 'df', 'all'...
%         'sm', 'df', 'all'};
%     labelArray = [row1; row2];
%     tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
%     ax = gca;
%     set(ax, 'box','off');
%     ax.LineWidth = 2;
%     ax.XTickLabel = tickLabels;
%     ax.XAxis.FontSize = 8;
%     ax.YAxis.FontSize = 20;
%     title(['Subject ' num2str(i)],'fontsize',20);
% end
% sgtitle('Two previous pattern matching','fontsize',30);
%
% %%
% figure();
% for i=1:(num_sub)
%     s=subplot(2,7,i);
%     cnt = 0;
%     for cs = 1:3
%         for cs_cs = 1:3
%             cnt = cnt + 1;
%             bar(cnt,squeeze(pattern_two_previous{i}(cs,cs_cs,1) - belief_two_previous{i}(cs,cs_cs,1)),'FaceColor',[0.75, 0.75, 0],'EdgeColor','k','LineWidth',0.75);
%             hold on;
%             errorbar(cnt,squeeze(pattern_two_previous{i}(cs,cs_cs,1) - belief_two_previous{i}(cs,cs_cs,1)),sqrt(squeeze(pattern_two_previous{i}(cs,cs_cs,2)^2 + belief_two_previous{i}(cs,cs_cs,2))^2),'-ko','LineWidth',2,'linestyle','none');
%             hold on;
%             yline(0.0,'k','LineWidth',2);
%             if i==1 || i==6
%                 ylabel('Probability','fontsize',30)
%             end
%         end
%     end
%     xlim([0.25 9.75]);
%     ylim([-1.0 1.0]);
%     xticks([1 2 3 4 5 6 7 8 9]);
%     row1 = {'2sc', '2sc', '2sc',...
%         '3sc', '3sc', '3sc'...
%         'all', 'all', 'all'};
%     row2 = {'sm', 'df', 'all',...
%         'sm', 'df', 'all'...
%         'sm', 'df', 'all'};
%     labelArray = [row1; row2];
%     tickLabels = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
%     ax = gca;
%     set(ax, 'box','off');
%     ax.LineWidth=2;
%     ax.XTickLabel = tickLabels;
%     ax.XAxis.FontSize = 8;
%     ax.YAxis.FontSize = 20;
%     title(['Subject ' num2str(i)],'fontsize',20);
% end
% sgtitle('Two previous difference between belief and pattern matching','fontsize',30);