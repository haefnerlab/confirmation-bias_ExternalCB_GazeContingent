% This is a script to display results of the analysis for the external CB project
% initialize parameters and hyperparameters
clear all; close all; clc;
boots = 50;
hi = 1.0;
low = 0.0;
points = 3;
hpr1 = 0.0;%logspace(-1, 5, 7);
hpr2 = logspace(-1, 5, 7);
test_thresh = 0;
frames = 4;
peri_cases = 1;
subjects = {'GCV3-subject01';'GCV3-subject02';'GCV3-subject03';'GCV3-subject07';'GCV3-subject08';'GCV3-subject11';'GCV3-subject13'};
[num_sub,~] = size(subjects);

not_wt_belief_extCB_opp = zeros(num_sub,points);
not_wt_belief_err_opp = zeros(num_sub,points);
not_wt_lobelief_err_opp = zeros(num_sub,points);
not_wt_hibelief_err_opp = zeros(num_sub,points);
not_wt_probability_chose_in_favor_opp = zeros(num_sub,points);


% pass on the parameters and subject ids to do respective analysis for each
% case of num_peri(number of elements in the periphery of the experiment: 2)
for sub=1:num_sub
    [params_boot, sobl, abbl, best_hprs,...
        data, accuracy, not_wt_belief_opp_cases, ...
        not_wt_trials_opp] = AllAnalysisGazeContingent_two_periphery_regeneration(subjects{sub}, 3, boots, points+1, hpr1, hpr2);
    [~,~,~,~,~,~, opp_c,opp_tr] = AllAnalysisGazeContingent_two_periphery_regeneration(subjects{sub}, 3, boots, 2, hpr1, hpr2);
    opp_cases(sub,:) = opp_c;
    opp_trials(sub) = opp_tr;
    sub_accuracy(sub) = accuracy;
    
    %assign PK weight related values from the returned results
    sub_opp_sacc(sub) = not_wt_trials_opp;
    trials_per_num_peri(sub) = (sum(data.num_peri(1:data.current_trial)==2));
    num_frames = 1 + frames * 2;
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
    
    not_wt_belief_mid_opp(sub,:) = not_wt_belief_opp_cases(1,:);
    not_wt_lo_err_opp(sub,:) = not_wt_belief_opp_cases(3,:);
    not_wt_hi_err_opp(sub,:) = not_wt_belief_opp_cases(4,:);
    not_wt_prob_chose_in_favor_opp(sub,:) = not_wt_belief_opp_cases(2,:);
    
end
%%
for sub=1:num_sub
    %     [belief_opp_cases, trials_opp,belief_opp_cases_fixed, trials_opp_fixed,not_wt_belief_opp_cases_fixed, not_wt_trials_opp_fixed, num_sessions, num_trials_per_session]
    [not_wt_divided_belief_opp_cases, not_wt_divided_trials_opp, params, num_sess, num_tr_pr_sess] = Data_Divided_AllAnalysisGazeContingent_two_periphery_regenerate(subjects{sub}, 3, boots, points+1,hpr1, hpr2);%squeeze(params_boot_fixed(sub,:,:))
    [opp_c_divided_temp, ~,~, ~, ~] = Data_Divided_AllAnalysisGazeContingent_two_periphery_regenerate(subjects{sub}, 3, boots, 2, hpr1, hpr2);
    num_session(sub) = num_sess;
    params_per_session{sub} = params;
    opp_cases_divided{sub} = opp_c_divided_temp;
    num_trials_per_session{sub} = num_tr_pr_sess;
    num_opp_sacc_per_session{sub} = not_wt_divided_trials_opp;
    not_wt_divided_belief_mid_opp{sub} = squeeze(not_wt_divided_belief_opp_cases(:,1,:));
    not_wt_divided_lo_err_opp{sub} = squeeze(not_wt_divided_belief_opp_cases(:,3,:));
    not_wt_divided_hi_err_opp{sub} = squeeze(not_wt_divided_belief_opp_cases(:,4,:));
    not_wt_divided_prob_chose_in_favor_opp{sub} = squeeze(not_wt_divided_belief_opp_cases(:,2,:));
end
%%
figure(1)
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

% this figure displays the results where the accumulated
% belief is compared to weighted signal in periphery
% figure(2)
% cnt = 1;
for i=1:(num_sub)
    disp(i)
    s=subplot(2,num_sub,num_sub+i);
    errorbar(squeeze(not_wt_belief_mid_opp(i,:)),squeeze(not_wt_prob_chose_in_favor_opp(i,:)),squeeze(not_wt_lo_err_opp(i,:)),squeeze(not_wt_hi_err_opp(i,:)),'-bo','LineWidth',2);
    hold on;
    yline(0.5,'k','LineWidth',2);
    xlabel('Accumulated evidence','Fontsize',15)
    ylabel('Prob chose in favor','Fontsize',15)
    %     title(['Ext CB for 2 in periphery for subject ' num2str(i)])
    %     xlim([0 350]);
    ylim([0 1.5]);
    ax = gca;
    ax.XAxis.FontSize = 12;
    ax.YAxis.FontSize = 12;
    text(s.XLim(1),s.YLim(2)-0.15,['Opp Saccade Num: ' num2str(floor(mean(sub_opp_sacc(i))))],'Fontsize',10);
end

%%
% this figure plots the accuracy of the subjects
figure(6)
subplot(1,2,1)
hold on;
bar(1:num_sub,sub_accuracy./trials_per_num_peri);
hold on;
yline(0.5,'k');
xlabel('Subject Number','Fontsize',20)
ylabel('Percent Correct','Fontsize',20)
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
title(['Accuracy under fixed noise and ratio of 0.7'],'Fontsize',18)


% this figure plots the accuracy of the subjects
% figure(7)
subplot(1,2,2)
hold on;
bar(1:num_sub,squeeze(opp_cases(:,2)));
hold on;
errorbar(1:num_sub,squeeze(opp_cases(:,2)),squeeze(opp_cases(:,3)),squeeze(opp_cases(:,4)),'ok')
xlabel('Subject Number','Fontsize',20)
ylabel('Prob chose in favor','Fontsize',20)
yline(0.5,'k');
ax = gca;
ax.XAxis.FontSize = 20;
ax.YAxis.FontSize = 20;
title(['Collapsed Prob chose in favor of accumulated evidence'],'Fontsize',18)

%%
% color = ['r','g','c','m','y'];
figure()
for i=1:(num_sub)
    disp(i)
    subplot(1,num_sub,i)
    hold on;
    for j=1:num_session(i)
        legend_id{j} = ['ses=' num2str(j) ',' 'sac=' num2str(num_opp_sacc_per_session{i}(j))];
        errorbar(squeeze(not_wt_divided_belief_mid_opp{i}(j,:)),squeeze(not_wt_divided_prob_chose_in_favor_opp{i}(j,:)),squeeze(not_wt_divided_lo_err_opp{i}(j,:)),squeeze(not_wt_divided_hi_err_opp{i}(j,:)),'-o','LineWidth',2);
        hold on;
    end
    hold on;
    legend_id{num_session(i)+1} = ['all sessions'];
    errorbar(squeeze(not_wt_belief_mid_opp(i,:)),squeeze(not_wt_prob_chose_in_favor_opp(i,:)),squeeze(not_wt_lo_err_opp(i,:)),squeeze(not_wt_hi_err_opp(i,:)),'-ko','LineWidth',2);
    hold on;
    legend(legend_id,'AutoUpdate','off');
    yline(0.5,'--k','LineWidth',2);
    ylim([0 1.5]);
    
    xlabel('Accumulated evidence','Fontsize',20)
    ylabel('Prob chose in favor','Fontsize',20)
end
%%
figure()
for i=1:(num_sub)
    disp(i)
    s(i)=subplot(2,num_sub,i)%num_sub+i)
    errorbar(1:num_session(i),squeeze(opp_cases_divided{i}(:,2)),squeeze(opp_cases_divided{i}(:,3)),squeeze(opp_cases_divided{i}(:,4)),'-o','LineWidth',2);
    %     hold on;
    %     errorbar(num_session(i)+1,squeeze(opp_cases(i,2)),squeeze(opp_cases(i,3)),squeeze(opp_cases(i,4)),'or')
    xlabel('Session Number')
    ylabel('Prob chose in favor')
    yline(0.5,'--k');
    ylim([0.3 1.0]);
    axis('tight');
    
end
linkaxes(s,'y');
for i=1:(num_sub)
    disp(i)
    s1(i)=subplot(2,num_sub,num_sub+i)
    plot(1:num_session(i),num_opp_sacc_per_session{i},'-ok','LineWidth',2);
    xlabel('Session Number')
    ylabel('Number of Opp categories in periphery')
    yline(10,'--k');
end
linkaxes(s1,'y');
%%

num_frames = 1 + frames * 2;
for i=1:(num_sub)
    figure()
    disp(i)
    subplot(2,num_sub,i)
    hold on;
    temporal_kernel_sess = [];
    norm_temporal_kernel_sess = [];
    lo_temporal_kernel_sess = [];
    hi_temporal_kernel_sess = [];
    lo_bias_sess = [];
    hi_bias_sess = [];
    for j=1:num_session(i)
        s(j)=subplot(1,num_session(i)+1,j)
        temporal_kernel_sess(j,:) = prctile(params_per_session{i}(j,:, 1:num_frames), 50);
        %         norm_temporal_kernel_sess(j,:) = temporal_kernel(j,:)/mean(temporal_kernel(i,:));
        lo_temporal_kernel_sess(j,:) = prctile(params_per_session{i}(j,:, 1:num_frames), 50) - prctile(params_per_session{i}(j,:, 1:num_frames), 16);
        hi_temporal_kernel_sess(j,:) = prctile(params_per_session{i}(j,:, 1:num_frames), 84) - prctile(params_per_session{i}(j,:, 1:num_frames), 50);
        lo_bias_sess(j) = prctile(params_per_session{i}(j,:, end-1), 50) - prctile(params_per_session{i}(j,:, end-1), 16);
        hi_bias_sess(j) = prctile(params_per_session{i}(j,:, end-1), 84) - prctile(params_per_session{i}(j,:, end-1), 50);
        
        errorbar(1:frames,squeeze(temporal_kernel_sess(j,1:frames)),squeeze(lo_temporal_kernel_sess(j,1:frames)),squeeze(hi_temporal_kernel_sess(j,1:frames)),'-ok','LineWidth',2)
        hold on;
        errorbar((1:frames-1),squeeze(temporal_kernel_sess(j,(frames+1):(end - pr ))),squeeze(lo_temporal_kernel_sess(j,(frames+1):(end - pr))),squeeze(hi_temporal_kernel_sess(j,(frames+1):(end - pr))),'or','LineStyle','none')
        hold on;
        errorbar([frames frames]',squeeze(temporal_kernel_sess(j,(end - pr+1 ):end)),squeeze(lo_temporal_kernel_sess(j,(end - pr+1):end)),squeeze(hi_temporal_kernel_sess(j,(end - pr+1):end)),'or','LineStyle','none')
        xlabel('Frames','Fontsize',20);
        ylabel('Weights','Fontsize',20);
        hold on;
        title(['Session:' num2str(j)])
        yline(0.0,'--k','LineWidth',2);
        hold on;
        %         errorbar([frames+1]',squeeze(prctile(params_per_session{i}(j,:, end-1), 50)),lo_bias_sess(j),hi_bias_sess(j),'ob','LineStyle','none')
    end
    hold on;
    s(num_session(i)+1)=subplot(1,num_session(i)+1,num_session(i)+1)
    errorbar(1:frames,squeeze(temporal_kernel(i,1:frames)),squeeze(lo_temporal_kernel(i,1:frames)),squeeze(hi_temporal_kernel(i,1:frames)),'-ok','LineWidth',2)
    hold on;
    errorbar((1:frames-1),squeeze(temporal_kernel(i,(frames+1):(end - pr ))),squeeze(lo_temporal_kernel(i,(frames+1):(end - pr))),squeeze(hi_temporal_kernel(i,(frames+1):(end - pr))),'or','LineStyle','none')
    hold on;
    errorbar([frames frames]',squeeze(temporal_kernel(i,(end - pr+1 ):end)),squeeze(lo_temporal_kernel(i,(end - pr+1):end)),squeeze(hi_temporal_kernel(i,(end - pr+1):end)),'or','LineStyle','none')
    hold on;
    %     errorbar([frames+1]',bias(i),bias_err(i,1),bias_err(i,2),'ob','LineStyle','none')
    yline(0.0,'--k','LineWidth',2);
    xlabel('Frames','Fontsize',20);
    ylabel('Weights','Fontsize',20);
    title(['Combined for all sessions'])
    linkaxes(s);
end

