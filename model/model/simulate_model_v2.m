function [chosen_locs,not_chosen_locs,final_choice]=simulate_model_v2(frame_signals,params,scale_normalize,is_sampling)

frame_signals=frame_signals/scale_normalize;

num_frames = 0.5*(size(frame_signals,2)-1);
num_trials = size(frame_signals,1);


num_samples_total_entropy=params(1);
num_samples_noise_entropy=params(2);
sigma_fovea = params(3);
sigma_peri = params(4);
p_match = params(5);
lapse_rate_choice=params(6);
lapse_bias_choice=params(7);
lapse_rate_saccade=params(8);
lapse_bias_saccade=params(9);

npts_quad = 5;
[pts_leg, wts_leg] = GaussLegendre(npts_quad);

pts_leg=norminv(0.5*pts_leg(:)+0.5,zeros(npts_quad,1),sigma_peri);
I_peri=cartprod(pts_leg,pts_leg,pts_leg,pts_leg,pts_leg,pts_leg,pts_leg,pts_leg);
wts_all=cartprod(wts_leg,wts_leg,wts_leg,wts_leg,wts_leg,wts_leg,wts_leg,wts_leg);
num_input=size(wts_all);
frame_signals=repmat(frame_signals,num_input,1);

rls_rt=(dec2bin(0:2^4-1)-'0')+1;
for i=1:size(rls_rt,1)
    mu1=pr_rt_il_rl(p_match,sigma_peri,sigma_fovea,[frame_signals(:,1),frame_signals(:,2)+I_peri(:,2),frame_signals(:,3)+I_peri(:,3)],[1,0,0]);
    
end


end

function pr=pr_rt_il_rl(p_match,sigma_peri,sigma_fovea,Il,rl)
Il(:,rl==1)=Il(:,rl==1)/(sigma_fovea.^2);
Il(:,rl==0)=Il(:,rl==0)/(sigma_peri.^2);
pr=sigmoid(exp(sum(log(p_match/(1-p_match))+log(sigmoid(2*Il-log(p_match/(1-p_match))))-log(sigmoid(2*Il+log(p_match/(1-p_match)))))));
end