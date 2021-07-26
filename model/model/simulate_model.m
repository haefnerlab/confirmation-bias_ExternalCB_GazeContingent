function [chosen_locs,not_chosen_locs,final_choice,log_odds_accumulated]=simulate_model(frame_signals,params,scale_normalize,is_sampling,orientation_std_exp)

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




l1_dim = 2;
l2_dim = 3;


prob_ct = zeros(num_trials,2,num_frames);
prob_ct(:,1,1) = 0.5;
prob_ct(:,2,1) = 0.5;
BAS = zeros(num_trials,num_frames-1,2);
table_ct_cl1_cl2_given_d = zeros(num_trials,2,2,2);
comp1 = zeros(num_trials,num_frames-1,2);
comp2 = zeros(num_trials,num_frames-1,2);
chosen_locs = zeros(num_trials,num_frames-1);
not_chosen_locs = zeros(num_trials,num_frames-1);
log_odds_accumulated = zeros(num_trials,num_frames+1);


final_choice = zeros(num_trials,1);
I_first = zeros(1,num_trials);
I_others = zeros(num_trials,2,num_frames);
I_peri = zeros(num_trials,2,num_frames);



for tr=1:num_trials
    if mod(tr,1000)==0
        disp(tr);
    end
    
    
    
    I_first(tr) = normrnd(frame_signals(tr,1),sigma_fovea0);
    I_others(tr,:,:) = normrnd([frame_signals(tr,2:2:end);frame_signals(tr,3:2:end)], sigma_fovea0);
    I_peri(tr,:,:) = normrnd([frame_signals(tr,2:2:end);frame_signals(tr,3:2:end)], sigma_peri0);
    
    
    log_odds_accumulated(tr,1) = logodds_model(I_first(tr),sigma_fovea.^2,p_match);
    
    
    for i=1:num_frames-1
        
        prob_ct(tr,1,i+1) = 0.5;
        prob_ct(tr,2,i+1) = 0.5;
        
        table_ct_cl1_cl2_given_d(tr,:,:,:) = compute_table(I_peri(tr,1,i),I_peri(tr,2,i),sigma_peri,squeeze(prob_ct(tr,:,i+1)),p_match,log_odds_accumulated(tr,i));
        
        sd1=randi(1e8);
        
        
        [BAS(tr,i,1),comp1(tr,i,1),comp2(tr,i,1)] = compute_score(squeeze(table_ct_cl1_cl2_given_d(tr,:,:,:)),l1_dim,num_samples_total_entropy,num_samples_noise_entropy,is_sampling,sd1);
        [BAS(tr,i,2),comp1(tr,i,2),comp2(tr,i,2)] = compute_score(squeeze(table_ct_cl1_cl2_given_d(tr,:,:,:)),l2_dim,num_samples_total_entropy,num_samples_noise_entropy,is_sampling,sd1);
        
        
        
        log_odds_l1=logodds_model(I_others(tr,1,i),sigma_fovea.^2,p_match);
        log_odds_l2=logodds_model(I_others(tr,2,i),sigma_fovea.^2,p_match);
        log_odds_l1_peri=logodds_model(I_peri(tr,1,i),sigma_peri.^2,p_match);
        log_odds_l2_peri=logodds_model(I_peri(tr,2,i),sigma_peri.^2,p_match);
        
        
        
        if rand()<lapse_rate_saccade
            if rand()<lapse_bias_saccade
                chosen_locs(tr,i) = 1;
                not_chosen_locs(tr,i) = 2;
                log_odds_accumulated(tr,i+1) = log_odds_accumulated(tr,i) + log_odds_l1 + log_odds_l2_peri + log_odds_l1_peri;
            else
                chosen_locs(tr,i) = 2;
                not_chosen_locs(tr,i) = 1;
                log_odds_accumulated(tr,i+1) = log_odds_accumulated(tr,i) + log_odds_l2 + log_odds_l1_peri + log_odds_l2_peri;
            end
            
        else
            
            if BAS(tr,i,1)>=BAS(tr,i,2)
                
                chosen_locs(tr,i) = 1;
                not_chosen_locs(tr,i) = 2;
                log_odds_accumulated(tr,i+1) = log_odds_accumulated(tr,i) + log_odds_l1 + log_odds_l2_peri + log_odds_l1_peri;
                
            else
                chosen_locs(tr,i) = 2;
                not_chosen_locs(tr,i) = 1;
                log_odds_accumulated(tr,i+1) = log_odds_accumulated(tr,i) + log_odds_l2 + log_odds_l1_peri + log_odds_l2_peri;
                
            end
        end
        
    end
    

    
    log_odds_l1_peri=logodds_model(I_peri(tr,1,end),sigma_peri.^2,p_match);
    log_odds_l2_peri=logodds_model(I_peri(tr,2,end),sigma_peri.^2,p_match);
    
    if isempty(i)
        i=0;
    end
    log_odds_accumulated(tr,end) = log_odds_accumulated(tr,i+1) + log_odds_l1_peri + log_odds_l2_peri;
    
    if rand()<lapse_rate_choice
        if rand()<lapse_bias_choice
            final_choice(tr) = 1;
        else
            final_choice(tr) = -1;
        end
    else
        if log_odds_accumulated(tr,end)>=0
            final_choice(tr) = -1;
        else
            final_choice(tr) = 1;
        end
    end
end
end
