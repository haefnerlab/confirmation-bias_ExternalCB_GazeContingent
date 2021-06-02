load('synthetic_data.mat')
sig_peri=logspace(-1,log10(5),14);
nss=[1:10,25,50,75,100];
orientation_std_exp=0.11;
for i=1:numel(nss)
    for j=1:numel(sig_peri)
        clc
        i
        j
        params=[100,nss(i),0.1,sig_peri(j),0.7,0,0.5,0,0.5];
        [chosen_locs,not_chosen_locs,final_choice]=simulate_model(frame_signals,params,scale_normalize,is_sampling,orientation_std_exp);
        hpr1 = logspace(-3, 3, 5);
        hpr2 = logspace(-3, 3, 5);
        frame_signals_ord0=frame_signals;
        for tr=1:num_trials
            frame_signals_ord0(tr,:)=[frame_signals(tr,1),frame_signals(tr,[1:2:2*(frames-1)]+chosen_locs(tr,:)),...
                frame_signals(tr,[1:2:2*(frames-1)]+not_chosen_locs(tr,:)),...
                frame_signals(tr,size(frame_signals,2)-1:size(frame_signals,2))];
        end
        frame_signals_ord=sign(frame_signals_ord0);
        trials = size(final_choice, 1);
        all_frames = num_images;
%         [best_hprs, ~] = xValidatePK_with_lapse(frame_signals_ord, final_choice(:)'==1, frames, hpr1, 0, hpr2, 0, 5);
%         [sig_wts, ~, ~, ~, ~, ~] = PsychophysicalKernelwithlapse(frame_signals_ord, final_choice==1, frames, best_hprs(1), 0, best_hprs(3), 0);
        [sig_wts, ~, ~, ~, ~, ~] = PsychophysicalKernelwithlapse(frame_signals_ord, final_choice(:)'==1, frames, 0, 0, 0, 0);
        % [sig_wts, ~, ~, ~, ~, ~] = PsychophysicalKernelwithlapse(frame_signals_ord, final_choice==1, frames, 0, 0, 0, 0);
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
        
        sig_wts1{i,j}=sig_wts;
        prop_confirmatory_saccades(i,j)=prop_confirmatory_saccade;
        accuracys(i,j)=mean(final_choice==category_trial);
    end
end