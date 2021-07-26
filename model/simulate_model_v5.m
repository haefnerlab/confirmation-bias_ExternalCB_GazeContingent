function [chosen_locs,not_chosen_locs,final_choice,log_odds_accumulated_agg]=simulate_model_v5(frame_signals,params,scale_normalize,is_sampling,orientation_std_exp)

frame_signals=frame_signals/scale_normalize;




num_frames = 0.5*(size(frame_signals,2)-1);
num_trials = size(frame_signals,1);
num_images = size(frame_signals,2);

num_samples_total_entropy=params(1);
num_samples_noise_entropy=params(2);

sigma_fovea = sqrt(orientation_std_exp.^2+params(3).^2);
sigma_fovea0 = params(3);

sigma_peri = sqrt(orientation_std_exp.^2+params(4).^2);
sigma_peri0 = params(4);

p_match = params(5);

im_fov=frame_signals+normrnd(0,sigma_fovea0,num_trials,num_images);
im_peri=frame_signals+normrnd(0,sigma_peri0,num_trials,num_images);




log_odds_accumulated = zeros(num_trials,num_frames+1);
log_odds_accumulated1 = zeros(num_trials,num_frames-1);


I_first = im_fov(:,1);
I_others = im_fov(:,2:end);
I_peri = im_peri(:,2:end);




log_odds_accumulated(:,1) = logodds_model(I_first,sigma_fovea.^2,p_match);


c1=zeros(num_trials,num_frames-1);
nc1=zeros(num_trials,num_frames-1);
gamma=(sigma_peri.^2)./((sigma_peri.^2)+(sigma_fovea.^2));
for i=1:num_frames-1
    
    sd1=randi(1e8);
    
    b1= compute_bas_v3(I_peri(:,2*i-1),sigmoid(-1*log_odds_accumulated(:,i)),p_match,sigma_peri^2,[num_samples_total_entropy,num_samples_noise_entropy],sd1);
    b2= compute_bas_v3(I_peri(:,2*i),sigmoid(-1*log_odds_accumulated(:,i)),p_match,sigma_peri^2,[num_samples_total_entropy,num_samples_noise_entropy],sd1);
%     log_odds_l1=logodds_model(I_others(:,2*i-1),sigma_fovea.^2,p_match);
%     log_odds_l2=logodds_model(I_others(:,2*i),sigma_fovea.^2,p_match);

    log_odds_l1=logodds_model(I_others(:,2*i-1)*gamma+I_peri(:,2*i-1)*(1-gamma),(sigma_fovea.^2)*gamma,p_match);
    log_odds_l2=logodds_model(I_others(:,2*i)*gamma+I_peri(:,2*i)*(1-gamma),(sigma_fovea.^2)*gamma,p_match);



    log_odds_l1_peri=logodds_model(I_peri(:,2*i-1),sigma_peri.^2,p_match);
    log_odds_l2_peri=logodds_model(I_peri(:,2*i),sigma_peri.^2,p_match);
    
    
    log_odds_accumulated1(:,i)=log_odds_accumulated(:,i);
    
    
    
    
    c1(b1>b2,i)=1;
    nc1(b1>b2,i)=2;
    log_odds_accumulated(b1>b2,i+1) = log_odds_accumulated(b1>b2,i) + log_odds_l1(b1>b2) + log_odds_l2_peri(b1>b2);
    
    c1(b1<b2,i)=2;
    nc1(b1<b2,i)=1;
    log_odds_accumulated(b1<b2,i+1) = log_odds_accumulated(b1<b2,i) + log_odds_l2(b1<b2) + log_odds_l2_peri(b1<b2);
    
    r1=rand(num_trials,1)<0.5;
    
    c1(r1==0& (b1==b2),i)=1;
    nc1(r1==0& (b1==b2),i)=2;
    log_odds_accumulated(r1==0 & (b1==b2),i+1) = log_odds_accumulated(r1==0& (b1==b2),i) + log_odds_l1(r1==0& (b1==b2)) + log_odds_l2_peri(r1==0& (b1==b2));
    
    
    c1(r1==1& (b1==b2),i)=2;
    nc1(r1==1& (b1==b2),i)=1;
    log_odds_accumulated(r1==1& (b1==b2),i+1) = log_odds_accumulated(r1==1& (b1==b2),i) + log_odds_l2(r1==1& (b1==b2)) + log_odds_l1_peri(r1==1& (b1==b2));
    
    
    
    
end
chosen_locs =c1;
not_chosen_locs = nc1;


log_odds_l1_peri=logodds_model(I_peri(:,end-1),sigma_peri.^2,p_match);
log_odds_l2_peri=logodds_model(I_peri(:,end),sigma_peri.^2,p_match);

if isempty(i)
    i=0;
end
log_odds_accumulated(:,end) = log_odds_accumulated(:,i+1) + log_odds_l1_peri + log_odds_l2_peri;
log_odds_accumulated_agg=log_odds_accumulated;
log_odds_accumulated_agg1=log_odds_accumulated1;


final_choice=1*ones(num_trials,1);
final_choice(log_odds_accumulated(:,end)>=0)=-1;

end
% end
