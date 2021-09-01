num_frames = 5;
num_trials = 1000;
p_match = 0.6;
c_t = sign(rand(num_trials, 1) - .5); %-ones(num_trials,1);

l1_dim_id = 2;
l2_dim_id = 3;
prob_ct = zeros(num_trials,2,num_frames);
prob_ct(:,1,1) = 0.5;
prob_ct(:,2,1) = 0.5;
BAS = zeros(num_trials,num_frames-1,2);

comp1 = zeros(num_trials,num_frames-1,2);
comp2 = zeros(num_trials,num_frames-1,2);
chosen_loc = zeros(num_trials,num_frames-1);
c_peri = zeros(num_trials,2,num_frames-1);
c_f1 = zeros(1,num_trials);


nt1=100;
log_odds_accumulated1 = zeros(nt1,num_trials,num_frames);
log_odds_frame_chosen1 = zeros(nt1,num_trials,num_frames);
log_odds_frame_not_chosen1 = zeros(nt1,num_trials,num_frames-1);
log_odds_frame_chosen_peri1 = zeros(nt1,num_trials,num_frames-1);
log_odds_frame_not_chosen_peri1 = zeros(nt1,num_trials,num_frames-1);


log_odds_accumulated = zeros(num_trials,num_frames);
log_odds_frame_chosen = zeros(num_trials,num_frames);
log_odds_frame_not_chosen = zeros(num_trials,num_frames-1);
log_odds_frame_chosen_peri = zeros(num_trials,num_frames-1);
log_odds_frame_not_chosen_peri = zeros(num_trials,num_frames-1);



choice = zeros(1,num_trials);
I_first = zeros(1,num_trials);
I_others = zeros(num_trials,2,num_frames-1);
I_peri = zeros(num_trials,2,num_frames-1);



ang_stim=deg2rad(45);
kappa_fix=1;
kappa_spatial=2;
[~,kappa_peri]=circ_vmpar(circ_vmrnd(0,kappa_fix,10000)+circ_vmrnd(0,kappa_spatial,10000));



nbins_I=5;
bins_I=linspace(-180,180,2+nbins_I);
bins_I(1)=[];
bins_I(end)=[];
bins_I=deg2rad(bins_I);


table_ct_Il1_Il2_given_d = zeros(num_trials,2,nbins_I,nbins_I);
for tr=1:num_trials
    if mod(tr,1000)==0
        
        disp(tr);
    end
    temp_peri = (rand(2,num_frames-1)<=p_match);%ones(2,num_frames-1);%
    c_peri(tr,:,:) = (2*temp_peri-1) * c_t(tr);
    
    c_f1(tr) = (2*(rand <= p_match)-1) * c_t(tr);
    ll=1;
    
    I0 = find_closest_binvalue(bins_I,(c_f1(tr)==-1)*circ_vmrnd(-1*ang_stim,kappa_fix,1)+(c_f1(tr)==1)*circ_vmrnd(ang_stim,kappa_fix,1));
    
    
    
    I_first(tr) = I0;
    I_others(tr,:,:) = find_closest_binvalue(bins_I,(squeeze(c_peri(tr,:,:))==-1).*circ_vmrnd(-ang_stim,kappa_fix,[size(c_peri,2),size(c_peri,3)]) + (squeeze(c_peri(tr,:,:))==1).*circ_vmrnd(ang_stim,kappa_fix,[size(c_peri,2),size(c_peri,3)]));
    I_peri(tr,:,:) = find_closest_binvalue(bins_I,squeeze(I_others(tr,:,:))+circ_vmrnd(0,kappa_spatial,[size(c_peri,2),size(c_peri,3)]));
    
    
    
    
    log_odds_accumulated1(ll,tr,1) = log(circ_vmpdf(I_first(tr),-ang_stim,kappa_fix)*p_match + circ_vmpdf(I_first(tr),ang_stim,kappa_fix)*(1-p_match)) ...
                                   - log(circ_vmpdf(I_first(tr),ang_stim,kappa_fix)*p_match + circ_vmpdf(I_first(tr),-ang_stim,kappa_fix)*(1-p_match));
    
    log_odds_frame_chosen1(ll,tr,1) = log_odds_accumulated1(ll,tr,1);
    for i=1:num_frames-1
        
        prob_ct(tr,1,i+1) = 0.5;
        prob_ct(tr,2,i+1) = 0.5;
        
        table_ct_Il1_Il2_given_d(tr,:,:,:) = compute_table_vm(I_peri(tr,1,i),I_peri(tr,2,i),kappa_fix,kappa_peri,prob_ct,p_match,sum(log_odds_accumulated1(ll,tr,:),3),nbins_I,bins_I,ang_stim);
        
        
        % Compute BAS scores (mutual informationb between c_t and category
        % at each location)
        [BAS(tr,i,1),comp1(tr,i,1),comp2(tr,i,1)] = compute_score_vm_v2(squeeze(table_ct_Il1_Il2_given_d(tr,:,:,:)),l1_dim_id,num_samples,is_sampling);
        [BAS(tr,i,2),comp1(tr,i,2),comp2(tr,i,2)] = compute_score_vm_v2(squeeze(table_ct_Il1_Il2_given_d(tr,:,:,:)),l2_dim_id,num_samples,is_sampling);
        %         fprintf('BAS(1) = %.2e\tBAS(2) = %.2e\n', BAS(tr,i,1), BAS(tr,i,2));
        
        % Compute log odds P(C_T = -1)/P(C_T = +1) for each location
        
        log_odds_l1 = log(circ_vmpdf(I_others(tr,1,i),-ang_stim,kappa_fix)*p_match+circ_vmpdf(I_others(tr,1,i),ang_stim,kappa_fix)*(1-p_match)) - log(circ_vmpdf(I_others(tr,1,i),ang_stim,kappa_fix)*p_match+circ_vmpdf(I_others(tr,1,i),-ang_stim,kappa_fix)*(1-p_match));
        log_odds_l2 = log(circ_vmpdf(I_others(tr,2,i),-ang_stim,kappa_fix)*p_match+circ_vmpdf(I_others(tr,2,i),ang_stim,kappa_fix)*(1-p_match)) - log(circ_vmpdf(I_others(tr,2,i),ang_stim,kappa_fix)*p_match+circ_vmpdf(I_others(tr,2,i),-ang_stim,kappa_fix)*(1-p_match));
        log_odds_l1_peri = log(circ_vmpdf(I_peri(tr,1,i),-ang_stim,kappa_peri)*p_match+circ_vmpdf(I_peri(tr,1,i),ang_stim,kappa_peri)*(1-p_match)) - log(circ_vmpdf(I_peri(tr,1,i),ang_stim,kappa_peri)*p_match+circ_vmpdf(I_peri(tr,1,i),-ang_stim,kappa_peri)*(1-p_match));
        log_odds_l2_peri = log(circ_vmpdf(I_peri(tr,2,i),-ang_stim,kappa_peri)*p_match+circ_vmpdf(I_peri(tr,2,i),ang_stim,kappa_peri)*(1-p_match)) - log(circ_vmpdf(I_peri(tr,2,i),ang_stim,kappa_peri)*p_match+circ_vmpdf(I_peri(tr,2,i),-ang_stim,kappa_peri)*(1-p_match));
        
        if BAS(tr,i,1)>=BAS(tr,i,2)
            chosen_loc(tr,i) = 1;
            I0 = I_others(tr,1,i);
            log_odds_accumulated1(ll,tr,i+1) = log_odds_accumulated1(ll,tr,i) + log_odds_l1 + log_odds_l2_peri;
            log_odds_frame_chosen1(ll,tr,i+1) = log_odds_l1;
            log_odds_frame_not_chosen1(ll,tr,i) = log_odds_l2;
            log_odds_frame_chosen_peri1(ll,tr,i) = log_odds_l1_peri;
            log_odds_frame_not_chosen_peri1(ll,tr,i) = log_odds_l2_peri;
        else
            chosen_loc(tr,i) = 2;
            I0 = I_others(tr,2,i);
            log_odds_accumulated1(ll,tr,i+1) = log_odds_accumulated1(ll,tr,i) + log_odds_l2 + log_odds_l1_peri;
            log_odds_frame_chosen1(ll,tr,i+1) = log_odds_l2;
            log_odds_frame_not_chosen1(ll,tr,i) = log_odds_l1;
            log_odds_frame_chosen_peri1(ll,tr,i) = log_odds_l2_peri;
            log_odds_frame_not_chosen_peri1(ll,tr,i) = log_odds_l1_peri;
        end
        
    end
    
    %     end
    log_odds_accumulated(tr,:)=mean(log_odds_accumulated1(:,tr,:),1);
    log_odds_frame_chosen(tr,:)=mean(log_odds_frame_chosen1(:,tr,:),1);
    log_odds_frame_not_chosen(tr,:)=mean(log_odds_frame_not_chosen1(:,tr,:),1);
    log_odds_frame_chosen_peri(tr,:)=mean(log_odds_frame_chosen_peri1(:,tr,:),1);
    log_odds_frame_not_chosen_peri(tr,:)=mean(log_odds_frame_not_chosen_peri1(:,tr,:),1);
    if log_odds_accumulated1(tr)>=0
        choice(tr) = -1;
    else
        choice(tr) = 1;
    end
end
function y = sigmoid(x)
y = (1 + exp(-x)).^-1;
end
