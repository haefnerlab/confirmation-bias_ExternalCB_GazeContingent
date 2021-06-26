% load('synthetic_data.mat')

scale_normalize=120;

orientation_std_exp=0.11;
rng('default')

a=[Inf];
b=[0.5];
% c=[1];
%
for i1=1:numel(a)
    for i2=1:numel(b)
%         clc
        i1
        i2
        
        %         for i3=1:numel(c)
        %             clc
        %             i1
        %             i2
        %             i3
        p_match=0.7;
        
        is_sampling=[1,1];
        params=[Inf,a(i1),0.1,b(i2),p_match,0,0.5,0,0.5];
        
        %%
        
        num_peri = 2;
        frames = 4;
        num_images = num_peri * frames +1;
        num_trials=1e4;
        
        
%         category_trial=sign(binornd(1,0.5,num_trials,1)-0.5);
        category_trial=sign(rand(num_trials,1)-0.5);
        frame_categories=repmat(category_trial,1,num_images).*sign(binornd(1,p_match,num_trials,num_images)-0.5);
        frame_signals=frame_categories+normrnd(0,orientation_std_exp,size(frame_categories));
        frame_signals=frame_signals*scale_normalize;
        
        
        %%
        % load('tmp');
        
        
        
        [chosen_locs,not_chosen_locs,final_choice,lo]=simulate_model_v5(frame_signals,params,scale_normalize,is_sampling,orientation_std_exp);
        
        
        %%
        num_reps=1;
        nn=1;
        
        hpr1 = logspace(-1, 5, 7);
        hpr2 = logspace(-1, 5, 7);
        
        
        
        for tr=1:num_trials
            frame_signals_ord0(tr,:)=[frame_signals(tr,1),frame_signals(tr,[1:2:2*(frames-1)]+chosen_locs(tr,:)),...
                frame_signals(tr,[1:2:2*(frames-1)]+not_chosen_locs(tr,:)),...
                frame_signals(tr,end-1:end)];
            
        end
        
        %             acc_bel=lo(:,1:(frames-1));
        %                 acc_bel=acc_bel(:);
        
        frame_signals_ord=sign(frame_signals_ord0);
        
        trials = size(final_choice, 1);
        all_frames = num_images;
        
        %                         [best_hprs, ~] = xValidatePK_with_lapse(frame_signals_ord, final_choice(:)'==1, frames, hpr1, 0, hpr2, 0, 10);
        %                         [sig_wts, ~, ~, ~, ~, ~] = PsychophysicalKernelwithlapse(frame_signals_ord, final_choice==1, frames, best_hprs(1), 0, best_hprs(3), 0);
        
        
%         [sig_wts, ~, ~, ~, ~, ~] = PsychophysicalKernelwithlapse(frame_signals_ord, final_choice(:)'==1, frames, 0, 0, 0, 0);
        [sig_wts, ~, ~, ~, ~, ~] = PsychophysicalKernelwithlapse(frame_signals_ord, final_choice(:)'==1, frames, 0.2955, 0, 0.4437, 0);
        
        
        
        % [sig_wts, ~, ~, ~, ~, ~] = PsychophysicalKernelwithlapse(frame_signals_ord, final_choice(:)'==1, frames, 0.29552, 0, 0.001, 0);
        weighted_signal_all = (sig_wts(1:all_frames))'.* frame_signals_ord;
        weighted_signal_chosen = weighted_signal_all(:,1:frames);
        weighted_signal_notchosen = weighted_signal_all(:,frames+1:(end-num_peri));
        
        accumulated_belief=weighted_signal_chosen(:,1);
        for i=1:(frames-2)
            accumulated_belief=[accumulated_belief,accumulated_belief(:,i)+weighted_signal_chosen(:,i+1)+weighted_signal_notchosen(:,i)];
        end
        
        % accumulated_belief=[weighted_signal_chosen(:,1),weighted_signal_chosen(:,1)+weighted_signal_chosen(:,2)+weighted_signal_notchosen(:,1)...
        %     ,weighted_signal_chosen(:,1)+weighted_signal_chosen(:,2)+weighted_signal_notchosen(:,1)+weighted_signal_chosen(:,3)+weighted_signal_notchosen(:,2)];
        sgn_signal_chosen = sign(frame_signals_ord(:,2:frames));
        sgn_signal_notchosen = sign(frame_signals_ord(:,frames+1:(end-num_peri)));
        all_sigs=[sgn_signal_chosen(:),sgn_signal_notchosen(:)];
        accumulated_belief=accumulated_belief(:);
        lo=lo(:);
        ids_del=prod(all_sigs,2)==1;
        all_sigs(ids_del,:)=[];
        accumulated_belief(ids_del)=[];
        lo(ids_del)=[];
        chosen_signal_sign=all_sigs(:,1);
        if frames>1
            prop_confirmatory_saccade(nn)=mean(chosen_signal_sign==sign(accumulated_belief));
            
        end
        
        accuracy(nn)=mean(final_choice==category_trial);
        wts_agg(:,nn)=sig_wts(1:end);
        
        
        %%
        
        parfor nn=2:num_reps
            clc
            nn
            
            category_trial=sign(binornd(1,0.5,num_trials,1)-0.5);
            frame_categories=repmat(category_trial,1,num_images).*sign(binornd(1,p_match,num_trials,num_images)-0.5);
            frame_signals=frame_categories+normrnd(0,orientation_std_exp,size(frame_categories));
            frame_signals=frame_signals*scale_normalize;
            
            
            
             [chosen_locs,not_chosen_locs,final_choice,lo]=simulate_model_v4(frame_signals,params,scale_normalize,is_sampling,orientation_std_exp);
            
            
            
            
            frame_signals_ord0=frame_signals;
            for tr=1:num_trials
                frame_signals_ord0(tr,:)=[frame_signals(tr,1),frame_signals(tr,[1:2:2*(frames-1)]+chosen_locs(tr,:)),...
                    frame_signals(tr,[1:2:2*(frames-1)]+not_chosen_locs(tr,:)),...
                    frame_signals(tr,end-1:end)];
                
            end
            
            frame_signals_ord=sign(frame_signals_ord0);
            trials = size(final_choice, 1);
            all_frames = num_images;
            
            
            
            [sig_wts, ~, ~, ~, ~, ~] = PsychophysicalKernelwithlapse(frame_signals_ord, final_choice(:)'==1, frames, 0, 0, 0, 0);
            %     [sig_wts, ~, ~, ~, ~, ~] = PsychophysicalKernelwithlapse(frame_signals_ord, final_choice(:)'==1, frames, 0.29552, 0, 0.001, 0);
            weighted_signal_all = (sig_wts(1:all_frames))'.* frame_signals_ord;
            weighted_signal_chosen = weighted_signal_all(:,1:frames);
            weighted_signal_notchosen = weighted_signal_all(:,(frames+1):(end-num_peri));
            
            accumulated_belief=weighted_signal_chosen(:,1);
            for i=1:(frames-2)
                accumulated_belief=[accumulated_belief,accumulated_belief(:,i)+weighted_signal_chosen(:,i+1)+weighted_signal_notchosen(:,i)];
            end
            
            % accumulated_belief=[weighted_signal_chosen(:,1),weighted_signal_chosen(:,1)+weighted_signal_chosen(:,2)+weighted_signal_notchosen(:,1)...
            %     ,weighted_signal_chosen(:,1)+weighted_signal_chosen(:,2)+weighted_signal_notchosen(:,1)+weighted_signal_chosen(:,3)+weighted_signal_notchosen(:,2)];
            
            
            sgn_signal_chosen = sign(frame_signals_ord(:,2:frames));
            sgn_signal_notchosen = sign(frame_signals_ord(:,frames+1:(end-num_peri)));
            all_sigs=[sgn_signal_chosen(:),sgn_signal_notchosen(:)];
            accumulated_belief=accumulated_belief(:);
            ids_del=prod(all_sigs,2)==1;
            all_sigs(ids_del,:)=[];
            accumulated_belief(ids_del)=[];
            chosen_signal_sign=all_sigs(:,1);
            if frames>1
                prop_confirmatory_saccade(nn)=mean(chosen_signal_sign==sign(accumulated_belief));
            end
            accuracy(nn)=mean(final_choice==category_trial);
            
            wts_agg(:,nn)=sig_wts(1:end);
            
            
            
        end
        wts_all{i1,i2}=wts_agg;
        acs{i1,i2}=accuracy;
        prs{i1,i2}=prop_confirmatory_saccade;
        
        %         end
    end
end

% ac1(jj)=mean(accuracy);
% if frames>1
%  pr1(jj)=mean(prop_confirmatory_saccade);
% end
% w1(:,jj)=mean(wts_agg,2);
