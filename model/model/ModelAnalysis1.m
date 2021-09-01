%% Data generation

num_peri = 2;
frames = 4;
num_images = num_peri * frames +1;
num_trials=300;
p_match=0.7;
orientation_stim=45;


hpr1 = logspace(-3, 3, 10);
hpr2 = logspace(-3, 3, 5);
boot_n = 500;

GaborData = newGaborData();

category_trial=sign(binornd(1,0.5,num_trials,1)-0.5);
frame_categories=repmat(category_trial,1,num_images).*sign(binornd(1,p_match,num_trials,num_images)-0.5);
frame_signals=zeros(size(frame_categories));
for tr=1:num_trials
    disp(tr);
    im = bpg.genImages(length(frame_categories(tr,:)), GaborData.stim_size, GaborData.stim_sp_freq_cpp,...
        GaborData.stim_std_sp_freq_cpp, frame_categories(tr,:)*orientation_stim, GaborData.noise(1), GaborData.annulus);
    
    im = uint8(im * GaborData.contrast(1) + 127);
    im = min(im, 255);
    im = max(im, 0);
    %Left category is 45 degree and right is -45. keeping it like this for
    %convention
    frame_signals(tr,:) = bpg.getSignal(double(im) - 127, GaborData.left_category, max(GaborData.noise(1), .04)) - ...
        bpg.getSignal(double(im) - 127, GaborData.right_category, max(GaborData.noise(1), .04));
end


clc


disp('Signals Generated');
save('synthetic_data2');
% %% Model simulation
% % load('synthetic_data');
% 
% % subid=3;
% % [frame_signals_ord_subj,frame_signals,final_choice_subj,category_trial]=load_data_subj(subid);
% % hpr1 = 0.0;
% % hpr2 = logspace(-1, 5, 7);
% % frames = 4;
% % num_peri=2;
% % p_match=0.7;
% % num_images = num_peri * frames +1;
% % num_trials=size(frame_signals,1);
% 
% scale_normalize=120;
% is_sampling=[1,1];
% 
% % num_samples_total_entropy=params(1);
% % num_samples_noise_entropy=params(2);
% % sigma_fovea = params(3);
% % sigma_peri = params(4);
% % p_match = params(5);
% % lapse_rate_choice=params(6);
% % lapse_bias_choice=params(7);
% % lapse_rate_saccade=params(8);
% % lapse_bias_saccade=params(9);
% 
% params=[100,100,0.1,0.3,p_match,0,0.5,0,0.5];
% [chosen_locs,not_chosen_locs,final_choice]=simulate_model(frame_signals,params,scale_normalize,is_sampling);
% disp('Model Simulated');
% 
% 
% 
% for tr=1:num_trials
%     frame_signals_ord0(tr,:)=[frame_signals(tr,1),frame_signals(tr,[1:2:2*(frames-1)]+chosen_locs(tr,:)),...
%         frame_signals(tr,[1:2:2*(frames-1)]+not_chosen_locs(tr,:)),...
%         frame_signals(tr,end-1:end)];
%     
% end
% 
% 
% 
% 
% %% Model predictions
% 
% accuracy=mean(final_choice==category_trial);
% 
% 
% 
% 
% [best_hprs, ~] = xValidatePK_with_lapse(frame_signals_ord, final_choice(:)'==1, frames, hpr1, 0, hpr2, 0, 10);
% disp('best params found')
% 
% 
% trials = size(final_choice, 1);
% all_frames = num_images;
% 
% [sig_wts, ~, ~, ~, ~, ~] = PsychophysicalKernelwithlapse(frame_signals_ord, final_choice==1, frames, best_hprs(1), 0, best_hprs(3), 0);
% 
% % signals weighted by PK weights based on final categorical choice
% weighted_signal_all = (sig_wts(1:all_frames))'.* frame_signals_ord;
% 
% weighted_signal_chosen = weighted_signal_all(:,1:frames);
% weighted_signal_notchosen = weighted_signal_all(:,frames+1:(end-num_peri));
% accumulated_belief=[weighted_signal_chosen(:,1),weighted_signal_chosen(:,1)+weighted_signal_chosen(:,2)+weighted_signal_notchosen(:,1)...
%     ,weighted_signal_chosen(:,1)+weighted_signal_chosen(:,2)+weighted_signal_notchosen(:,1)+weighted_signal_chosen(:,3)+weighted_signal_notchosen(:,2)];
% % weighted_signal_notchosen = reshape(weighted_signal_notchosen,trials,num_peri-1,frames-1);
% 
% 
% % ideal signals not weighted at all
% sgn_signal_chosen = sign(frame_signals_ord(:,2:frames));
% sgn_signal_notchosen = sign(frame_signals_ord(:,frames+1:(end-num_peri)));
% % signal_notchosen = reshape(signal_notchosen,trials,num_peri-1,frames-1);
% 
% 
% all_sigs=[sgn_signal_chosen(:),sgn_signal_notchosen(:)];
% accumulated_belief=accumulated_belief(:);
% ids_del=prod(all_sigs,2)==1;
% all_sigs(ids_del,:)=[];
% accumulated_belief(ids_del)=[];
% chosen_signal_sign=all_sigs(:,1);
% 
% prop_confirmatory_saccade=mean(chosen_signal_sign==sign(accumulated_belief));
% 
% %% Subject predictions
% 
% 
% accuracy_subj=mean(final_choice_subj==category_trial);
% 
% 
% 
% [best_hprs, ~] = xValidatePK_with_lapse(frame_signals_ord_subj, final_choice_subj(:)'==1, frames, hpr1, 0, hpr2, 0, 10);
% disp('best params found')
% 
% 
% trials = size(final_choice_subj, 1);
% all_frames = num_images;
% 
% [sig_wts, ~, ~, ~, ~, ~] = PsychophysicalKernelwithlapse(frame_signals_ord_subj, final_choice_subj==1, frames, best_hprs(1), 0, best_hprs(3), 0);
% 
% % signals weighted by PK weights based on final categorical choice
% weighted_signal_all = (sig_wts(1:all_frames))'.* frame_signals_ord_subj;
% 
% weighted_signal_chosen = weighted_signal_all(:,1:frames);
% weighted_signal_notchosen = weighted_signal_all(:,frames+1:(end-num_peri));
% accumulated_belief=[weighted_signal_chosen(:,1),weighted_signal_chosen(:,1)+weighted_signal_chosen(:,2)+weighted_signal_notchosen(:,1)...
%     ,weighted_signal_chosen(:,1)+weighted_signal_chosen(:,2)+weighted_signal_notchosen(:,1)+weighted_signal_chosen(:,3)+weighted_signal_notchosen(:,2)];
% % weighted_signal_notchosen = reshape(weighted_signal_notchosen,trials,num_peri-1,frames-1);
% 
% 
% % ideal signals not weighted at all
% sgn_signal_chosen = sign(frame_signals_ord_subj(:,2:frames));
% sgn_signal_notchosen = sign(frame_signals_ord_subj(:,frames+1:(end-num_peri)));
% % signal_notchosen = reshape(signal_notchosen,trials,num_peri-1,frames-1);
% 
% 
% all_sigs=[sgn_signal_chosen(:),sgn_signal_notchosen(:)];
% accumulated_belief=accumulated_belief(:);
% ids_del=prod(all_sigs,2)==1;
% all_sigs(ids_del,:)=[];
% accumulated_belief(ids_del)=[];
% chosen_signal_sign=all_sigs(:,1);
% 
% prop_confirmatory_saccade_subj=mean(chosen_signal_sign==sign(accumulated_belief));
% 
% 
% 
% %%
% 
% figure(1);
% subplot(2,2,1)
% hline(0.5,'b--');
% hold on
% bar(2*subid-1,accuracy_subj,0.5,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0,0,0],'LineWidth',1.5);
% bar(2*subid,accuracy,0.5,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0,0,0],'LineWidth',1.5);
% errorbar(2*subid-1,accuracy_subj,accuracy_subj.*(1-accuracy_subj)./sqrt(num_trials),'ko','linewidth',2)
% errorbar(2*subid,accuracy,accuracy.*(1-accuracy)./sqrt(num_trials),'ko','linewidth',2)
% xticks([1:2*subid]);
% xlab={};
% for i=1:subid
%     xlab=[xlab,'Data','Model'];
% end
% xticklabels(xlab);
% 
% ylim([0,1]);
% xlabel('Subject ID');
% ylabel('Accuracy');
% set(gca,'fontsize',20,'fontweight','bold')
% set(gca,'linewidth',4)
% box off
% 
% 
% 
% figure(1);
% subplot(2,2,2)
% hline(0.5,'b--');
% hold on
% bar(2*subid-1,prop_confirmatory_saccade_subj,0.5,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0,0,0],'LineWidth',1.5);
% bar(2*subid,prop_confirmatory_saccade,0.5,'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0,0,0],'LineWidth',1.5);
% errorbar(2*subid-1,prop_confirmatory_saccade_subj,prop_confirmatory_saccade_subj.*(1-prop_confirmatory_saccade_subj)./sqrt(num_trials),'ko','linewidth',2)
% errorbar(2*subid,prop_confirmatory_saccade,prop_confirmatory_saccade.*(1-prop_confirmatory_saccade)./sqrt(num_trials),'ko','linewidth',2)
% xticks([1:2*subid]);
% xticklabels(xlab);
% ylim([0,1]);
% xlabel('Subject ID');
% ylabel({'prop. confirmatory','saccade'});
% set(gca,'fontsize',20,'fontweight','bold')
% set(gca,'linewidth',4)
% box off
% 
% 
% set(gcf, 'color', 'white')