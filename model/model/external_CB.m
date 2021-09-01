clear;
colors = ['b','g','r','k','c','m','y'];


conditions{1} = [0 1 1; 0 100 100;0 1 2];
conditions{2} = [0 1 1; 0 1 5;0 1 2];

for jj=1:2
    
    ck = 1;
    for icond=1:size(conditions{jj}, 2)
        is_sampling = repmat(conditions{jj}(1, icond),1,2);
        num_samples_total_entropy = conditions{jj}(2, icond);
        num_samples_noise_entropy = conditions{jj}(3, icond);
        eye_movement_simulation
        num_elems = 5;
        mn_sig = min(log_odds_accumulated(:));
        mx_sig = max(log_odds_accumulated(:));
        ll = max(abs(mn_sig),abs(mx_sig));
        t1 = log_odds_accumulated(:,1:end-2);
        has_both_categories = squeeze(sum(c_peri(:,:,1:end-1), 2)) == 0;
        choice_sign = zeros(num_trials, num_frames-1);
        for tr=1:num_trials
            for fr=1:num_frames-1
                choice_sign(tr, fr) = c_peri(tr, chosen_loc(tr, fr), fr);
            end
        end
        has_both_categories = has_both_categories(:);
        t1 =t1(:);
        choice_sign=choice_sign(:);
        t1_when_both_choices = t1(has_both_categories);
        [num_vals, bin_edges, bin_idx] = histcounts(t1_when_both_choices, num_elems);
        bin_centers = (bin_edges(2:end) + bin_edges(1:end-1))/2;
        choice_sign_when_both = choice_sign(has_both_categories);
        
        for ibin=1:length(bin_centers)
            [prob_chose_neg(ibin), prob_chose_neg_CI(:, ibin)] = ...
                binofit(sum(choice_sign_when_both(bin_idx == ibin) == -1), num_vals(ibin));
            [prob_chose_pos(ibin), prob_chose_pos_CI(:, ibin)] = ...
                binofit(sum(choice_sign_when_both(bin_idx == ibin) == 1), num_vals(ibin));
        end
        
        if is_sampling
            displayname = sprintf('# Samples (component 1) = %d\n# Samples (component 2) = %d', num_samples_total_entropy,num_samples_noise_entropy);
        else
            displayname = 'Ideal';
        end
        
        figure(jj);
        %     subplot(2,2,jj)
        hold on;
        
        
        errorbar(bin_centers, prob_chose_neg, (prob_chose_neg - prob_chose_neg_CI(1, :)), (prob_chose_neg_CI(2,:) - prob_chose_neg), 'o-', 'Color', colors(ck), 'DisplayName', displayname, 'LineWidth', 2);
        hline(0.5,'k--');
        set(gca,'fontsize',20,'fontweight','bold');
        xlabel('Accumulated evidence for C_T','fontsize',20,'fontweight','bold');
        ylabel({'Saccade probability to','evidence in favor of category C_T'},'fontsize',20,'fontweight','bold');
%         title('Approximating entropy computation of both terms','fontsize',20,'fontweight','bold');
        %     title('Approximating entropy computation only in the conditional entropy','fontsize',20,'fontweight','bold');
        legend('Location', 'best');
        legend boxoff
        
        set(gca,'linewidth',2)
        ylim([0,1]);
        
        ck = ck+1;
        drawnow;
        
        summary{icond}.cat_trials=c_t>0;          %0 or 1
        summary{icond}.cat_first_frame=c_f1(:)*(45);
        summary{icond}.cat_frames=c_peri*(45);
        summary{icond}.chosen_loc=chosen_loc;
        summary{icond}.final_choice=choice(:);
        
    end
    set(gcf,'color','white')
    box off
end
