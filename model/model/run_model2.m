load('synthetic_data.mat')
params=[100,1,0.1,0.3,0.7,0,0.5,0,0.5];

[chosen_locs,not_chosen_locs,final_choice]=simulate_model_v4(frame_signals,params,scale_normalize,is_sampling,0);

%%
% load('tmp');
hpr1 = logspace(-3, 3, 10);
hpr2 = 0;%logspace(-3, 3, 5);
for tr=1:num_trials
    frame_signals_ord0(tr,:)=[frame_signals(tr,1),frame_signals(tr,1+chosen_locs(tr,1)),frame_signals(tr,3+chosen_locs(tr,2)),frame_signals(tr,5+chosen_locs(tr,3)),...
        frame_signals(tr,1+not_chosen_locs(tr,1)),frame_signals(tr,3+not_chosen_locs(tr,2)),frame_signals(tr,5+not_chosen_locs(tr,3)),...
        frame_signals(tr,end-1:end)];
    
end

frame_signals_ord=frame_signals_ord0;
trials = size(final_choice, 1);
all_frames = num_images;

% [best_hprs, ~] = xValidatePK_with_lapse(frame_signals_ord, final_choice(:)'==1, frames, hpr1, 0, hpr2, 0, 5);
% [sig_wts, ~, ~, ~, ~, ~] = PsychophysicalKernelwithlapse(frame_signals_ord, final_choice==1, frames, best_hprs(1), 0, best_hprs(3), 0);
[sig_wts, ~, ~, ~, ~, ~] = PsychophysicalKernelwithlapse(frame_signals_ord, final_choice==1, frames, 0, 0, 0, 0);
weighted_signal_all = (sig_wts(1:all_frames))'.* frame_signals_ord;
weighted_signal_chosen = weighted_signal_all(:,1:frames);
weighted_signal_notchosen = weighted_signal_all(:,frames+1:(end-num_peri));
accumulated_belief=[weighted_signal_chosen(:,1),weighted_signal_chosen(:,1)+weighted_signal_chosen(:,2)+weighted_signal_notchosen(:,1)...
    ,weighted_signal_chosen(:,1)+weighted_signal_chosen(:,2)+weighted_signal_notchosen(:,1)+weighted_signal_chosen(:,3)+weighted_signal_notchosen(:,2)];
sgn_signal_chosen = sign(frame_signals_ord(:,2:frames));
sgn_signal_notchosen = sign(frame_signals_ord(:,frames+1:(end-num_peri)));
all_sigs=[sgn_signal_chosen(:),sgn_signal_notchosen(:)];
accumulated_belief=accumulated_belief(:);
ids_del=prod(all_sigs,2)==1;
all_sigs(ids_del,:)=[];
accumulated_belief(ids_del)=[];
chosen_signal_sign=all_sigs(:,1);
prop_confirmatory_saccade=mean(chosen_signal_sign==sign(accumulated_belief));
accuracy=mean(final_choice==category_trial);