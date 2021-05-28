load('synthetic_data.mat')
nss=[1,3,100];
nboot=100;
hpr1 = logspace(-3, 3, 10);
hpr2 = logspace(-3, 3, 10);

for i=1:numel(nss)
    params=[100,nss(i),0.1,0.27,0.7,0,0.5,0,0.5];
    [chosen_locs,not_chosen_locs,final_choice]=simulate_model(frame_signals,params,scale_normalize,is_sampling);
    for tr=1:num_trials
        frame_signals_ord0(tr,:)=[frame_signals(tr,1),frame_signals(tr,1+chosen_locs(tr,1)),frame_signals(tr,3+chosen_locs(tr,2)),frame_signals(tr,5+chosen_locs(tr,3)),...
            frame_signals(tr,1+not_chosen_locs(tr,1)),frame_signals(tr,3+not_chosen_locs(tr,2)),frame_signals(tr,5+not_chosen_locs(tr,3)),...
            frame_signals(tr,end-1:end)];
    end
    frame_signals_ord=frame_signals_ord0;
    trials = size(final_choice, 1);
    all_frames = num_images;
    
    [best_hprs, ~] = xValidatePK_with_lapse(frame_signals_ord, final_choice(:)'==1, frames, hpr1, 0, hpr2, 0, 5);
    [sig_wts, ~, ~, ~, ~, ~] = PsychophysicalKernelwithlapse(frame_signals_ord, final_choice==1, frames, best_hprs(1), 0, best_hprs(3), 0);
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
    prop_confirmatory_saccade{i}(1)=mean(chosen_signal_sign==sign(accumulated_belief));
    accuracy{i}(1)=mean(final_choice==category_trial);
    ch_sn{i,1}=chosen_signal_sign(:);
    ac_bf{i,1}=accumulated_belief(:);
    wts{i}(:,1)=sig_wts(:);
    
    [q1,q2,q3]=histcounts(ac_bf{i,1},10);
    for k=1:numel(q1)
        xs{i}(k,1)=0.5*(q2(k)+q2(k+1));
        ys{i}(k,1)=mean(1==sign(ch_sn{i,1}(q3==k)));
        
    end
    
    frame_signals_ord0=frame_signals_ord;
    final_choice0=final_choice;
    for j=1:nboot
        clc
        i
        j
        ids=randi(num_trials,[num_trials,1]);
        frame_signals_ord=frame_signals_ord0(ids,:);
        final_choice=final_choice0(ids);
%         [best_hprs, ~] = xValidatePK_with_lapse(frame_signals_ord, final_choice(:)'==1, frames, hpr1, 0, hpr2, 0, 5);
        [sig_wts, ~, ~, ~, ~, ~] = PsychophysicalKernelwithlapse(frame_signals_ord, final_choice==1, frames, best_hprs(1), 0, best_hprs(3), 0);
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
        prop_confirmatory_saccade{i}(1+j)=mean(chosen_signal_sign==sign(accumulated_belief));
        accuracy{i}(1+j)=mean(final_choice==category_trial);
        ch_sn{i,1+j}=chosen_signal_sign(:);
        ac_bf{i,1+j}=accumulated_belief(:);
        wts{i}(:,1+j)=sig_wts(:);
        
        [q1,q2,q3]=histcounts(ac_bf{i,1+j},10);
        for k=1:numel(q1)
            xs{i}(k,j+1)=0.5*(q2(k)+q2(k+1));
            ys{i}(k,j+1)=mean(1==sign(ch_sn{i,1+j}(q3==k)));
            
        end
        
        
    end
    
end