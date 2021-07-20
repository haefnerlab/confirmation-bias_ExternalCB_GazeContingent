load('GazeContingentExternalCB_02-Jun-2021.mat');
load('belief_baseline.mat');
fx = @(mu,nt,baseline_bias) normcdf(abs(mu-baseline_bias)./sqrt((mu*(1-mu))/nt),'upper');
for sub=1:num_sub
    for pp=1:length(points)
    baseline_bias_sub{pp}(sub,:) = get_baseline_belief(squeeze(belief_mid_opp{pp}(sub,:)));
    end
    p_values_opp_cases(sub) = fx(prob_chose_in_favor_opp{1}(sub,:),sub_opp_sacc(sub),baseline_bias_sub{1}(sub));
    if (p_values_opp_cases(sub))<0.001
        stars{sub} = '***';
    elseif (p_values_opp_cases(sub))<0.01
        stars{sub} = '**';
    elseif (p_values_opp_cases(sub))<0.05
        stars{sub} = '*';
    else
        stars{sub} = ' ';
    end
end

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
bar(1:num_sub,squeeze(prob_chose_in_favor_opp{1})-baseline_bias_sub{1}, 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor','k','LineWidth',0.75);
bar(authors,squeeze(prob_chose_in_favor_opp{1}(authors))-baseline_bias_sub{1}(authors), 'FaceColor',[0.9290 0.6940 0.1250],'EdgeColor',[0.6350 0.0780 0.1840],'LineWidth',3);
hold on;
errorbar(1:num_sub,squeeze(prob_chose_in_favor_opp{1})-baseline_bias_sub{1},squeeze(lo_err_opp{1}),squeeze(hi_err_opp{1}),'ok','LineWidth', 1.5)
hold on;
text([1:num_sub]-0.3,(squeeze(prob_chose_in_favor_opp{1})-baseline_bias_sub{1})+sign(squeeze(prob_chose_in_favor_opp{1})-baseline_bias_sub{1})*0.065,stars,'FontSize',15,'FontWeight','bold');
xlabel('Subject Number','Fontsize',20)
ylabel({'Probability of saccading to ';'confirming evidence'},'Fontsize',20);
ax = gca;
ax.LineWidth=2;
ylim([-0.2 0.5]);
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
for pp=2:length(points)
    figure();
    set(gcf, 'renderer', 'painters');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperSize', [4 2]);
    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperPosition', [0 0 4 2]);
    set(gca, 'Position', get(gca, 'OuterPosition') - ...
    get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
    for i=1:(num_sub)
        s=subplot(2,sb_plt,i);
        errorbar(squeeze(belief_mid_opp{pp}(i,:)),squeeze(prob_chose_in_favor_opp{pp}(i,:)),squeeze(lo_err_opp{pp}(i,:)),squeeze(hi_err_opp{pp}(i,:)),'-bo','LineWidth',2);
        hold on;
        plot(squeeze(belief_mid_opp{pp}(i,:)),squeeze(baseline_bias_sub{pp}(i,:)),'-ro','LineWidth',2);
%         errorbar(squeeze(opp_raw_sig{pp}(i,:)),squeeze(opp_raw_prob{pp}(i,:)),squeeze(opp_raw_lo_err{pp}(i,:)),squeeze(opp_raw_hi_err{pp}(i,:)),'-ro','LineWidth',2);
        yline(0.5,'k','LineWidth',2);
        if i==1 || i==(sb_plt+1)
            ylabel('Probability chose in favor','fontsize',30);
        end
        if i==md_x || i==(sb_plt+md_x)
            xlabel('Accumulated evidence','fontsize',30);
        end
        yline(0.5,'k','LineWidth',2);
        xticks(round(squeeze(belief_mid_opp{pp}(i,:)),1));
        xlim([min(belief_mid_opp{pp}(i,:))-0.01 max(squeeze(belief_mid_opp{pp}(i,:)))+0.01]);
        ylim([0.3 1.0]);
%         text(belief_mid_opp{pp}(i,1)+0.1,0.2,['Opp. Sacc: ' num2str(floor(mean(sub_opp_sacc(i))))],'fontsize',15);
        ax = gca;
        set(ax, 'box','off');
        ax.LineWidth=2;
        ax.XAxis.FontSize = 20;
        ax.YAxis.FontSize = 20;
        title(['Subject ' num2str(i)],'fontsize',20);
    end
    sgtitle(['Bias w.r.t accumulated evidence for ' num2str(pp) ' bins'],'fontsize',30);
end
