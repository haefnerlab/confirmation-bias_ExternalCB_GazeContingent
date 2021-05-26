% This is a script to display results of the analysis for the external CB project
% initialize parameters and hyperparameters
clear all; close all; clc;
warning('off');
sanity_bins = 10;
boots = 50;
points = [1 2 3]; % number of points in accumulated belief vs bias plot
ps_points = 3;
hpr1 = logspace(-1, 5, 7);
hpr2 = logspace(-1, 5, 7);
frames = 4; % number of fixations
subjects = ...
    {'GCV3-subject01';'GCV3-subject02';...
    'GCV3-subject03';'GCV3-subject07';...
    'GCV3-subject08';'GCV3-subject11';...
    'GCV3-subject13';'GCV3-subject15';...
    'GCV3-subject16';'GCV3-subject17'};
[num_sub,~] = size(subjects);

for i=1:length(points)
    not_wt_belief_extCB_opp{i} = zeros(num_sub,points(i));
    not_wt_belief_err_opp{i} = zeros(num_sub,points(i));
    not_wt_lobelief_err_opp{i} = zeros(num_sub,points(i));
    not_wt_hibelief_err_opp{i} = zeros(num_sub,points(i));
    not_wt_probability_chose_in_favor_opp{i} = zeros(num_sub,points(i));
end

fx=@(mu,nt) normcdf(abs(mu-0.5)./sqrt((mu*(1-mu))/nt),'upper');

% pass on the parameters and subject ids to do respective analysis for each
% case of num_peri(number of elements in the periphery of the experiment: 2)
for sub=1:num_sub
    disp(['Starting analysis for subject ' num2str(sub) ' ...']);
    tic;
    [params_boot, sobl, abbl, best_hprs,...
        data, accuracy, not_wt_belief_opp_cases, ...
        not_wt_trials_opp,opp_r,...
        mn_sig(sub,:),mn_sig_norm(sub,:),mn_sig_ideal_foveal_only(sub,:),mn_sig_ideal_foveal_only_norm(sub,:),...
        mn_sig_ideal_all(sub,:),mn_sig_ideal_all_norm(sub,:),mn_perf(sub,:),mn_perf_ideal_foveal_only(sub,:),mn_perf_ideal_all(sub,:),...
        tr_seg(sub,:),err_perf(sub,:),tr_seg_ideal_foveal_only(sub,:),err_perf_ideal_foveal_only(sub,:),tr_seg_ideal_all(sub,:),err_perf_ideal_all(sub,:),...
        sanity_actual_data(sub,:,:), sanity_predicted_data(sub,:,:), sanity_mean_sig_shown(sub,:)]...
        = AllAnalysisGazeContingent_two_periphery_regeneration(subjects{sub}, 3, boots, points+1, hpr1, hpr2, ps_points, sanity_bins);
    
    
    %assign PK weight related values from the returned results
    sub_opp_sacc(sub) = not_wt_trials_opp{1};
    trials_per_num_peri(sub) = (sum(data.num_peri(1:data.current_trial)==2));
    num_frames = 1 + frames * 2;
%     alpha(sub,:) = [prctile(params_boot(:,end).^2,50) sqrt(var(params_boot(:,end).^2)/size(params_boot,1))];%[prctile(params_boot(:, end).^2, 50) std(params_boot(:, end).^2)];
    alpha(sub,:) = [prctile(1e-4+(1-1e-4)*sigmoid(params_boot(:,end)),50) sqrt(var(1e-4+(1-1e-4)*sigmoid(params_boot(:,end)))/size(params_boot,1))];%[prctile(params_boot(:, end).^2, 50) std(params_boot(:, end).^2)];
    bias(sub) = prctile(exp(params_boot(:, end-1)), 50);
    temporal_kernel(sub,:) = prctile(exp(params_boot(:, 1:num_frames)), 50);
    norm_temporal_kernel(sub,:) = temporal_kernel(sub,:)/mean(temporal_kernel(sub,:));
    lo_temporal_kernel(sub,:) = prctile(params_boot(:, 1:num_frames), 50) - prctile(params_boot(:, 1:num_frames), 16);
    hi_temporal_kernel(sub,:) = prctile(params_boot(:, 1:num_frames), 84) - prctile(params_boot(:, 1:num_frames), 50);
    beta(sub) = prctile(squeeze(abbl(:,2)),50);
    norm_all_linear(sub,:,:) = [sobl(:,1)/mean(temporal_kernel(sub,1:frames)) sobl(:,2)/mean(temporal_kernel(sub,1:frames)) sobl(:,3) sobl(:,4)];%sobl_norm;
    norm_slope_all(sub,:) = norm_all_linear(sub,:,1);
    norm_slope(sub) = prctile(squeeze(norm_all_linear(sub,:,1)),50);
    hprs_used(sub,:) = best_hprs;
   
    for pp=1:length(points)
        not_wt_belief_mid_opp{pp}(sub,:) = not_wt_belief_opp_cases{pp}(1,:);
        not_wt_lo_err_opp{pp}(sub,:) = not_wt_belief_opp_cases{pp}(3,:);
        not_wt_hi_err_opp{pp}(sub,:) = not_wt_belief_opp_cases{pp}(4,:);
        not_wt_prob_chose_in_favor_opp{pp}(sub,:) = not_wt_belief_opp_cases{pp}(2,:);
        opp_raw_sig{pp}(sub,:) = opp_r{pp}(1,:);
        opp_raw_prob{pp}(sub,:) = opp_r{pp}(2,:);
        opp_raw_lo_err{pp}(sub,:) = opp_r{pp}(3,:);
        opp_raw_hi_err{pp}(sub,:) = opp_r{pp}(4,:);
    end
    
    opp_trials(sub) = not_wt_trials_opp{1};
    sub_accuracy(sub) = sum(accuracy);
    
    p_values_opp_cases(sub) = fx(not_wt_prob_chose_in_favor_opp{1}(sub,:),opp_trials(sub));
    if (p_values_opp_cases(sub))<0.001
        stars{sub} = '***';
    elseif (p_values_opp_cases(sub))<0.01
        stars{sub} = '**';
    elseif (p_values_opp_cases(sub))<0.05
        stars{sub} = '*';
    else
        stars{sub} = ' ';
    end
    disp(['Completed analysis for subject ' num2str(sub) ' !!']);
    disp(['The learned regression weights: ' num2str(temporal_kernel(sub,:))]);
    disp(['The learned lapse: ' num2str(alpha(sub,1))]);
    toc;
    disp('-----------------------------------------------------------------------------------------------------');
    
    
    figure();
    subplot(1,2,1)
    pr=2;
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
    ax.LineWidth=2;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title(['Subject ' num2str(sub)],'fontsize',20);
    pause(1.5);
    
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
    ax.LineWidth=2;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
    title(['Subject ' num2str(sub)],'fontsize',20);
    pause(1.5);
end

%%
fig2 = figure();
set(fig2,'defaultLegendAutoUpdate','off');
for sub=1:(num_sub)
    s(sub)=subplot(2,5,sub);
    if sub==1 || sub==6
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
    if sub==3 || sub==8
        xlabel('Signal Displayed','fontsize',30);
    end
    if sub==1 || sub==6
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
%%
fig = figure();
set(fig,'defaultLegendAutoUpdate','off');
pr=2;
for i=1:(num_sub)
    s(i)=subplot(2,5,i);
    if i==1 || i==6
        LH(1)=errorbar(1:frames,squeeze(temporal_kernel(i,1:frames)),squeeze(lo_temporal_kernel(i,1:frames)),squeeze(hi_temporal_kernel(i,1:frames)),'-ob','LineWidth',2);
        L{1} = 'Foveated';
        hold on;
        LH(2)= errorbar((2:frames),squeeze(temporal_kernel(i,(frames+1):(end - pr ))),squeeze(lo_temporal_kernel(i,(frames+1):(end - pr))),squeeze(hi_temporal_kernel(i,(frames+1):(end - pr))),'or','LineStyle','none','LineWidth',2);
        L{2} = 'Not foveated';
        hold on;
    else
        errorbar(1:frames,squeeze(temporal_kernel(i,1:frames)),squeeze(lo_temporal_kernel(i,1:frames)),squeeze(hi_temporal_kernel(i,1:frames)),'-ob','LineWidth',2);
        hold on;
        errorbar((2:frames),squeeze(temporal_kernel(i,(frames+1):(end - pr ))),squeeze(lo_temporal_kernel(i,(frames+1):(end - pr))),squeeze(hi_temporal_kernel(i,(frames+1):(end - pr))),'or','LineStyle','none','LineWidth',2);
        hold on;
        errorbar([frames+1 frames+1]',squeeze(temporal_kernel(i,(end - pr+1 ):end)),squeeze(lo_temporal_kernel(i,(end - pr+1):end)),squeeze(hi_temporal_kernel(i,(end - pr+1):end)),'or','LineStyle','none','LineWidth',2);
    end
    if i==3 || i==8
        xlabel('Fixation Number','fontsize',30);
    end
    if i==1 || i==6
        ylabel('Unnormalized Weights','fontsize',30);
    end
    if i==1
        legend(LH,L, 'Fontsize',20, 'Box','off');
    end
    hold('on');
    xlim([1 frames + 1]);
    text(2,2,[num2str(trials_per_num_peri(i)) ' trials'],'fontsize',15);
    yline(0,'k','linewidth',2);
%     ylim([-0.2 0.3]);
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
    s=subplot(2,5,i);
    bar(1,squeeze(not_wt_prob_chose_in_favor_opp{1}(i,:)),'FaceColor','b','EdgeColor','k','LineWidth',0.75);
    hold on;
    errorbar(1,squeeze(not_wt_prob_chose_in_favor_opp{1}(i,:)),squeeze(not_wt_lo_err_opp{1}(i,:)),squeeze(not_wt_hi_err_opp{1}(i,:)),'-ko','LineWidth',2,'linestyle','none');
    hold on;
    yline(0.5,'k','LineWidth',2);
    if i==1 || i==6
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
        s=subplot(2,5,i);
        errorbar(squeeze(not_wt_belief_mid_opp{pp}(i,:)),squeeze(not_wt_prob_chose_in_favor_opp{pp}(i,:)),squeeze(not_wt_lo_err_opp{pp}(i,:)),squeeze(not_wt_hi_err_opp{pp}(i,:)),'-bo','LineWidth',2);
        hold on;
        %     errorbar(squeeze(opp_raw_sig{pp}(i,:)),squeeze(opp_raw_prob{pp}(i,:)),squeeze(opp_raw_lo_err{pp}(i,:)),squeeze(opp_raw_hi_err{pp}(i,:)),'-ro','LineWidth',2);
        yline(0.5,'k','LineWidth',2);
        if i==1 || i==6
            ylabel('Probability chose in favor','fontsize',30);
        end
        if i==3 || i==8
            xlabel('Accumulated evidence','fontsize',30);
        end
        yline(0.5,'k','LineWidth',2);
%         xlim([0 50]);
        ylim([0.2 1.0]);
        text(not_wt_belief_mid_opp{pp}(i,1),0.4,['Opp. Sacc: ' num2str(floor(mean(sub_opp_sacc(i))))],'fontsize',15);
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
subplot(2,3,1)
hold on;
bar(1:num_sub,sub_accuracy./trials_per_num_peri, 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',0.75);
hold on;
bar([2 3],sub_accuracy([2 3])./trials_per_num_peri([2 3]), 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor',[0.6350 0.0780 0.1840],'LineWidth',3.0);
hold on;
yline(0.5,'k','LineWidth',1.5);
xlabel('Subject Number','Fontsize',20);
ylabel('Percent Correct','Fontsize',20);
ax = gca;
ax.LineWidth=2;
ylim([0.4 1]);
xlim([0.5 10.5]);
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
ax.XTick = [1 2 3 4 5 6 7 8 9 10];
hold on;
tmp = sub_accuracy./trials_per_num_peri;
errorbar(1:num_sub,tmp,sqrt((tmp.*(1-tmp))./trials_per_num_peri),'ok','LineWidth', 1.5);
hold on;

subplot(2,3,2)
hold on;
bar(1:num_sub,squeeze(not_wt_prob_chose_in_favor_opp{1}), 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',0.75);
bar([2 3],squeeze(not_wt_prob_chose_in_favor_opp{1}([2 3])), 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor',[0.6350 0.0780 0.1840],'LineWidth',3);
hold on;
errorbar(1:num_sub,squeeze(not_wt_prob_chose_in_favor_opp{1}),squeeze(not_wt_lo_err_opp{1}),squeeze(not_wt_hi_err_opp{1}),'ok','LineWidth', 1.5)
hold on;
text([1:num_sub]-0.3,squeeze(not_wt_prob_chose_in_favor_opp{1})+0.065,stars,'FontSize',15,'FontWeight','bold');
xlabel('Subject Number','Fontsize',20)
ylabel({'Probability of saccading to ';'confirming evidence'},'Fontsize',20);
yline(0.5,'k','LineWidth',1.5);
ax = gca;
ax.LineWidth=2;
ylim([0.4 1]);
xlim([0.5 10.5]);
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
ax.XTick = [1 2 3 4 5 6 7 8 9 10];
hold on;

subplot(2,3,3)
for i=1:(num_sub)
    if i==2 || i==3
        plot(1:frames,squeeze(temporal_kernel(i,1:frames)),'--k');
    else
        plot(1:frames,squeeze(temporal_kernel(i,1:frames)),'k');
    end
    hold on;
end
LH(1) = plot(1:frames,mean(squeeze(temporal_kernel(:,1:frames)),1),'-ok','LineWidth',2);
L{1} = 'Foveated';
ylim([-0.05 0.05]);
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
    if i==2 || i==3
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
legend(LH,L, 'Fontsize',20, 'Box','off');

sgtitle('All subject all analysis','fontsize',30);
%%
figure();
k1 = 1;
k2 = 2;
k = [];
for sub=1:num_sub
    if sub==num_sub
        LH6(1) = bar(k1, squeeze(not_wt_prob_chose_in_favor_opp{1}(sub,:)),'FaceColor',[.5 0 .5],'EdgeColor','k','LineWidth',0.75);
        L6{1} = 'Weighted fovea + periphery';
    else
        bar(k1, squeeze(not_wt_prob_chose_in_favor_opp{1}(sub,:)),'FaceColor',[.5 0 .5],'EdgeColor','k','LineWidth',0.75);
    end
    hold on;
    if sub==num_sub
        LH6(2) = bar(k2, opp_raw_prob{1}(sub,:),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.75);
        L6{2} = 'Fovea only';
    else
        bar(k2, opp_raw_prob{1}(sub,:),'FaceColor',[0.8500 0.3250 0.0980],'EdgeColor','k','LineWidth',0.75);
    end
    hold on;
    errorbar([k1 k2],[squeeze(not_wt_prob_chose_in_favor_opp{1}(sub,:)) opp_raw_prob{1}(sub,:)],[squeeze(not_wt_lo_err_opp{1}(sub,:)) opp_raw_lo_err{1}(sub,:)],[squeeze(not_wt_hi_err_opp{1}(sub,:)) opp_raw_hi_err{1}(sub,:)],'ok','LineWidth',2,'linestyle','none');
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
title('Sanity Check: Left: Weighted fovea + periphery, Right: Fovea only','fontsize',30);

%%
figure();
subplot(1,2,1)
cnt = 1;
for i=1:(num_sub)
    if i==2||i==3
        plot([1 2],squeeze(not_wt_prob_chose_in_favor_opp{2}(i,:)),'-ro','LineWidth',1.0);
    else
        plot([1 2],squeeze(not_wt_prob_chose_in_favor_opp{2}(i,:)),'-bo','LineWidth',1.0);
    end
    hold on;
end
hold on;
scatter(ones(1,num_sub),squeeze(not_wt_prob_chose_in_favor_opp{2}(:,1)),100,'b','filled');
hold on;
scatter(ones(1,2),squeeze(not_wt_prob_chose_in_favor_opp{2}([2 3],1)),100,'r','filled');
hold on;
LH(1) = scatter(ones(1,num_sub)*2,squeeze(not_wt_prob_chose_in_favor_opp{2}(:,2)),100,'b','filled');
L{1} = 'Naive Subject';
hold on;
LH(2) = scatter(ones(1,2)*2,squeeze(not_wt_prob_chose_in_favor_opp{2}([2 3],2)),100,'r','filled');
L{2} = 'Author';
hold on;
plot([1 2],squeeze(mean(not_wt_prob_chose_in_favor_opp{2},1)),'-ko','LineWidth',5);
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
errorbar(squeeze(not_wt_prob_chose_in_favor_opp{2}(:,1)),squeeze(not_wt_prob_chose_in_favor_opp{2}(:,2)),squeeze(not_wt_lo_err_opp{2}(:,1)),squeeze(not_wt_hi_err_opp{2}(:,1)),'horizontal','LineWidth',2, 'LineStyle', 'none','Color','k');
hold on;
errorbar(squeeze(not_wt_prob_chose_in_favor_opp{2}(:,1)),squeeze(not_wt_prob_chose_in_favor_opp{2}(:,2)),squeeze(not_wt_lo_err_opp{2}(:,2)),squeeze(not_wt_hi_err_opp{2}(:,2)),'vertical', 'LineWidth',2,'LineStyle', 'none','Color','k');
hold on;
scatter(squeeze(not_wt_prob_chose_in_favor_opp{2}(:,1)),squeeze(not_wt_prob_chose_in_favor_opp{2}(:,2)),100,'b','filled');
hold on;
scatter(squeeze(not_wt_prob_chose_in_favor_opp{2}([2 3],1)),squeeze(not_wt_prob_chose_in_favor_opp{2}([2 3],2)),100,'r','filled');
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
hold on;

%%
chosen_sub = [5 7];
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
bar(1:length(chosen_sub),squeeze(squeeze(not_wt_prob_chose_in_favor_opp{1}(chosen_sub,:))), 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',0.75);
hold on;
errorbar(1:length(chosen_sub),squeeze(squeeze(not_wt_prob_chose_in_favor_opp{1}(chosen_sub,:))),squeeze(not_wt_lo_err_opp{1}(chosen_sub,:)),squeeze(not_wt_hi_err_opp{1}(chosen_sub,:)),'ok','LineWidth', 1.5);
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
ylim([-0.05 0.05]);
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
%%
figure();
for i=1:num_sub
    subplot(2,5,i)
    errorbar(mn_sig_ideal_all_norm(i,:),mn_perf_ideal_all(i,:),err_perf_ideal_all(i,:),'og','linewidth',2);
    hold on;
    errorbar(mn_sig_norm(i,:),mn_perf(i,:),err_perf(i,:),'ob','linewidth',2);
    hold on;
    plot(mn_sig_ideal_foveal_only_norm(i,:),mn_perf_ideal_foveal_only(i,:),'-or','linewidth',2);
    hold on;
    errorbar(mn_sig_ideal_foveal_only_norm(i,:),mn_perf_ideal_foveal_only(i,:),err_perf_ideal_foveal_only(i,:),'or','linewidth',2);
    hold on;
    plot(mn_sig_ideal_all_norm(i,:),mn_perf_ideal_all(i,:),'-og','linewidth',2);
    hold on;
    plot(mn_sig_norm(i,:),mn_perf(i,:),'-ob','linewidth',2);
    title(['Subject ' num2str(i)],'fontsize',20);
    if i==3 || i==8
        xlabel('Normalized Evidence Accumulated','fontsize',20);
    end
    if i==1 || i==6
        ylabel('Probability of correct answer','fontsize',20);
    end
    ylim([0.25 1.001]);
    xlim([0.0 1.001]);
    ax = gca;
    set(ax, 'box','off');
    ax.LineWidth=2;
    ax.XAxis.FontSize = 20;
    ax.YAxis.FontSize = 20;
end

sgtitle(['Green: Ideal all frames Red: Ideal fovea only Blue: Subject'], 'fontsize',30);