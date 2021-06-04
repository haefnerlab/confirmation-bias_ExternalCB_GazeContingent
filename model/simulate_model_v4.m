function [chosen_locs,not_chosen_locs,final_choice,log_odds_accumulated_agg1]=simulate_model_v4(frame_signals,params,scale_normalize,is_sampling,orientation_std_exp)

frame_signals=frame_signals/scale_normalize;




num_frames = 0.5*(size(frame_signals,2)-1);
num_trials = size(frame_signals,1);


num_samples_total_entropy=params(1);
num_samples_noise_entropy=params(2);

sigma_fovea = sqrt(orientation_std_exp.^2+params(3).^2);
sigma_fovea0 = params(3);

sigma_peri = sqrt(orientation_std_exp.^2+params(4).^2);
sigma_peri0 = params(4);

p_match = params(5);
lapse_rate_choice=params(6);
lapse_bias_choice=params(7);
lapse_rate_saccade=params(8);
lapse_bias_saccade=params(9);


chosen_locs = zeros(num_trials,num_frames-1);
not_chosen_locs = zeros(num_trials,num_frames-1);
log_odds_accumulated_agg=zeros(num_trials,num_frames+1);
log_odds_accumulated_agg1=zeros(num_trials,num_frames-1);

final_choice = zeros(num_trials,1);



parfor tr=1:num_trials
    
    log_odds_accumulated = zeros(1,num_frames+1);
    log_odds_accumulated1 = zeros(1,num_frames-1);
%     if mod(tr,1000)==0
%         disp(tr);
%     end
    
    I_first = normrnd(frame_signals(tr,1),sigma_fovea0);
    I_others = normrnd([frame_signals(tr,2:2:end);frame_signals(tr,3:2:end)], sigma_fovea0);
    I_peri = normrnd([frame_signals(tr,2:2:end);frame_signals(tr,3:2:end)], sigma_peri0);
    
    
    log_odds_accumulated(1) = logodds_model(I_first,sigma_fovea.^2,p_match);
    
    
    c1=zeros(1,num_frames-1);
    nc1=zeros(1,num_frames-1);
    for i=1:num_frames-1
        
        sd1=randi(1e8);
        
        b1= compute_bas(I_peri(1,i),sigmoid(-1*log_odds_accumulated(i)),p_match,sigma_peri^2,[num_samples_total_entropy,num_samples_noise_entropy],sd1);
        b2= compute_bas(I_peri(2,i),sigmoid(-1*log_odds_accumulated(i)),p_match,sigma_peri^2,[num_samples_total_entropy,num_samples_noise_entropy],sd1);
        
        log_odds_l1=logodds_model(I_others(1,i),sigma_fovea.^2,p_match);
        log_odds_l2=logodds_model(I_others(2,i),sigma_fovea.^2,p_match);
        log_odds_l1_peri=logodds_model(I_peri(1,i),sigma_peri.^2,p_match);
        log_odds_l2_peri=logodds_model(I_peri(2,i),sigma_peri.^2,p_match);
        
        
        log_odds_accumulated1(i)=log_odds_accumulated(i);%+ log_odds_l2_peri + log_odds_l1_peri;
        if rand()<lapse_rate_saccade
            if rand()<lapse_bias_saccade
                c1(i) = 1;
                nc1(i) = 2;
                log_odds_accumulated(i+1) = log_odds_accumulated(i) + log_odds_l1 + log_odds_l2_peri + log_odds_l1_peri;
            else
                c1(i) = 2;
                nc1(i) = 1;
                log_odds_accumulated(i+1) = log_odds_accumulated(i) + log_odds_l2 + log_odds_l1_peri + log_odds_l2_peri;
            end
            
        else
            
            if b1>=b2
                
                c1(i) = 1;
                nc1(i) = 2;
                log_odds_accumulated(i+1) = log_odds_accumulated(i) + log_odds_l1 + log_odds_l2_peri + log_odds_l1_peri;
                
            else
                c1(i) = 2;
                nc1(i) = 1;
                log_odds_accumulated(i+1) = log_odds_accumulated(i) + log_odds_l2 + log_odds_l1_peri + log_odds_l2_peri;
                
            end
        end
        
    end
    chosen_locs(tr,:) =c1;
    not_chosen_locs(tr,:) = nc1;
    
    
    log_odds_l1_peri=logodds_model(I_peri(1,end),sigma_peri.^2,p_match);
    log_odds_l2_peri=logodds_model(I_peri(2,end),sigma_peri.^2,p_match);
    
    if isempty(i)
        i=0;
    end
    log_odds_accumulated(end) = log_odds_accumulated(i+1) + log_odds_l1_peri + log_odds_l2_peri;
    log_odds_accumulated_agg(tr,:)=log_odds_accumulated;
    log_odds_accumulated_agg1(tr,:)=log_odds_accumulated1;
    if rand()<lapse_rate_choice
        if rand()<lapse_bias_choice
            final_choice(tr) = 1;
        else
            final_choice(tr) = -1;
        end
    else
        if log_odds_accumulated(end)>=0
            final_choice(tr) = -1;
        else
            final_choice(tr) = 1;
        end
    end
end
end
