num_frames = 4;
num_trials = 1000;
p_match = 0.7;
c_t = sign(rand(num_trials, 1) -0.5);
sigma0 = 0.01;
sigma1 = 0.01;
ecc_sig = 0.4;
total_sigma_periph = sqrt(sigma0^2+ecc_sig^2);
l1_dim = 2;
l2_dim = 3;
prob_ct = zeros(num_trials,2,num_frames);
prob_ct(:,1,1) = 0.5;
prob_ct(:,2,1) = 0.5;
BAS = zeros(num_trials,num_frames-1,2);
table_ct_cl1_cl2_given_d = zeros(num_trials,2,2,2);
comp1 = zeros(num_trials,num_frames-1,2);
comp2 = zeros(num_trials,num_frames-1,2);
chosen_loc = zeros(num_trials,num_frames-1);
c_peri = zeros(num_trials,2,num_frames);
c_f1 = zeros(1,num_trials);



log_odds_accumulated = zeros(num_trials,num_frames+1);
% log_odds_frame_chosen = zeros(num_trials,num_frames);
% log_odds_frame_not_chosen = zeros(num_trials,num_frames-1);
% log_odds_frame_chosen_peri = zeros(num_trials,num_frames-1);
% log_odds_frame_not_chosen_peri = zeros(num_trials,num_frames-1);






choice = zeros(1,num_trials);
I_first = zeros(1,num_trials);
I_others = zeros(num_trials,2,num_frames);
I_peri = zeros(num_trials,2,num_frames);



for tr=1:num_trials
    if mod(tr,1000)==0
        disp(tr);
    end
    temp_peri = (rand(2,num_frames)<=p_match);
    c_peri(tr,:,:) = (2*temp_peri-1) * c_t(tr);
    c_f1(tr) = (2*(rand <= p_match)-1) * c_t(tr);
    
    I0 = normrnd(1*c_f1(tr),sigma0);
    I_first(tr) = I0;
    I_others(tr,:,:) = normrnd(squeeze(1*c_peri(tr,:,:)), sigma0);
    I_peri(tr,:,:) = normrnd(squeeze(I_others(tr,:,:)), sqrt(sigma1.^2+ecc_sig.^2));
    I_others(tr,:,:) = normrnd(squeeze(I_others(tr,:,:)), sigma1);
    log_odds_accumulated(tr,1) = log(normpdf(I_first(tr),-1,sigma0)*p_match + normpdf(I_first(tr),1,sigma0)*(1-p_match)) ...
        - log(normpdf(I_first(tr),1,sigma0)*p_match + normpdf(I_first(tr),-1,sigma0)*(1-p_match));
    %     log_odds_frame_chosen(tr,1) = log_odds_accumulated(tr,1);
    for i=1:num_frames-1
        
        prob_ct(tr,1,i+1) = 0.5;
        prob_ct(tr,2,i+1) = 0.5;
        
        table_ct_cl1_cl2_given_d(tr,:,:,:) = compute_table(I_peri(tr,1,i),I_peri(tr,2,i),total_sigma_periph,squeeze(prob_ct(tr,:,i+1)),p_match,log_odds_accumulated(tr,i));
        
        sd1=randi(1e8);
        
        
        [BAS(tr,i,1),comp1(tr,i,1),comp2(tr,i,1)] = compute_score(squeeze(table_ct_cl1_cl2_given_d(tr,:,:,:,:)),l1_dim,num_samples_total_entropy,num_samples_noise_entropy,is_sampling,sd1);
        [BAS(tr,i,2),comp1(tr,i,2),comp2(tr,i,2)] = compute_score(squeeze(table_ct_cl1_cl2_given_d(tr,:,:,:,:)),l2_dim,num_samples_total_entropy,num_samples_noise_entropy,is_sampling,sd1);
        
        log_odds_l1 = log(normpdf(I_others(tr,1,i),-1,sigma0)*p_match + normpdf(I_others(tr,1,i),1,sigma0)*(1-p_match)) - log(normpdf(I_others(tr,1,i),-1,sigma0)*(1-p_match) + normpdf(I_others(tr,1,i),1,sigma0)*p_match);
        log_odds_l2 = log(normpdf(I_others(tr,2,i),-1,sigma0)*p_match + normpdf(I_others(tr,2,i),1,sigma0)*(1-p_match)) - log(normpdf(I_others(tr,2,i),-1,sigma0)*(1-p_match) + normpdf(I_others(tr,2,i),1,sigma0)*p_match);
        log_odds_l1_peri = log(normpdf(I_peri(tr,1,i),-1,total_sigma_periph)*p_match + normpdf(I_peri(tr,1,i),1,total_sigma_periph)*(1-p_match)) - log(normpdf(I_peri(tr,1,i),-1,total_sigma_periph)*(1-p_match) + normpdf(I_peri(tr,1,i),1,total_sigma_periph)*p_match);
        log_odds_l2_peri = log(normpdf(I_peri(tr,2,i),-1,total_sigma_periph)*p_match + normpdf(I_peri(tr,2,i),1,total_sigma_periph)*(1-p_match)) - log(normpdf(I_peri(tr,2,i),-1,total_sigma_periph)*(1-p_match) + normpdf(I_peri(tr,2,i),1,total_sigma_periph)*p_match);
        
        if BAS(tr,i,1)>=BAS(tr,i,2)
            chosen_loc(tr,i) = 1;
            
            log_odds_accumulated(tr,i+1) = log_odds_accumulated(tr,i) + log_odds_l1 + log_odds_l2_peri + log_odds_l1_peri;
            %             log_odds_frame_chosen(tr,i+1) = log_odds_l1;
            %             log_odds_frame_not_chosen(tr,i) = log_odds_l2;
            %             log_odds_frame_chosen_peri(tr,i) = log_odds_l1_peri;
            %             log_odds_frame_not_chosen_peri(tr,i) = log_odds_l2_peri;
        else
            chosen_loc(tr,i) = 2;
            
            log_odds_accumulated(tr,i+1) = log_odds_accumulated(tr,i) + log_odds_l2 + log_odds_l1_peri + log_odds_l2_peri;
            %             log_odds_frame_chosen(tr,i+1) = log_odds_l2;
            %             log_odds_frame_not_chosen(tr,i) = log_odds_l1;
            %             log_odds_frame_chosen_peri(tr,i) = log_odds_l2_peri;
            %             log_odds_frame_not_chosen_peri(tr,i) = log_odds_l1_peri;
        end
        
    end
    log_odds_l1_peri = log(normpdf(I_peri(tr,1,end),-1,total_sigma_periph)*p_match + normpdf(I_peri(tr,1,end),1,total_sigma_periph)*(1-p_match)) - log(normpdf(I_peri(tr,1,end),-1,total_sigma_periph)*(1-p_match) + normpdf(I_peri(tr,1,end),1,total_sigma_periph)*p_match);
    log_odds_l2_peri = log(normpdf(I_peri(tr,2,end),-1,total_sigma_periph)*p_match + normpdf(I_peri(tr,2,end),1,total_sigma_periph)*(1-p_match)) - log(normpdf(I_peri(tr,2,end),-1,total_sigma_periph)*(1-p_match) + normpdf(I_peri(tr,2,end),1,total_sigma_periph)*p_match);
    log_odds_accumulated(tr,end) = log_odds_accumulated(tr,i) + log_odds_l1_peri + log_odds_l2_peri;
    
    if log_odds_accumulated(tr,end)>=0
        choice(tr) = 0;
    else
        choice(tr) = 1;
    end
end

