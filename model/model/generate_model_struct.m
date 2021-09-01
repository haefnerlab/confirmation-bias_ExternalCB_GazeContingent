 function model = generate_model_struct(sub_struct,incomplete_dir_save)
model.params_names = {...
    'Beta (choice)',...
    'Beta (saccade)',...
    'Sensory noise (fovea)',...
    'Sensory noise (periphery)',...
    'Number of samples (total entropy)',...
    'Number of samples (noise entropy)',...
    'p_match',...
    'Lapse rate',...
    'Lapse bias'...
    };


model.id = id;


if exist('incomplete_dir_save','var')
    model.incomplete_dir_save=incomplete_dir_save;
end
% Parameter family:
% Beta(a,b)=1
% Gamma(k,theta)=2
% Normal(mu,sig)=3
% Lognormal(mu,sig)=4
% Beta-prime(mu,sig)=5





model.min_sig=1e-4;
model.min_lapse=1e-4;

model.params_family = [1,1,2,2,1,1,1,1,1];
model.params_additive_offset = [1e-4,1e-4,1e-4,1e-4,1,1,0.5,1e-4,1e-4];
model.params_multiplicative_offset = [(1-1e-4),(1-1e-4),1,1,99,99,0.5,(1-1e-4),(1-1e-4)];
model.params_hyperprior_params{1} = [1.25,1.25,0.5,0.5,1,1,2,1,1.25];
model.params_hyperprior_params{2} = [1.25,1.25,1,1,1,1,3,5,1.25];


model.exact_inference=0;
model.normprior_std = 1;
model.max_phi=6*model.normprior_std;
model.ns_inf=100;
model.eps_max=atand((0.5*60.96)/50); %Maximum eccentricity

model.num_params = length(model.params_family);
model.num_eff_params = model.num_params;

model.phi_theta = @(phi, model) phi_theta(phi, model);
model.theta_phi = @(theta, model) theta_phi(theta, model);

model.npts_quad = 20;
[model.pts_leg, model.wts_leg] = GaussLegendre(model.npts_quad);



model.design_matrix{1}.num_frames=sub_struct.num_frames;
model.design_matrix{1}.num_frames=sub_struct.num_frames;



model.set_default = [];
model.default_values = [];
[~, ids_sort] = sort(model.set_default, 'descend');
model.num_eff_params = model.num_params - numel(model.set_default);
model.set_default = model.set_default(ids_sort);
model.default_values = model.default_values(ids_sort);

model.params_family_eff = model.params_family;
model.params_additive_offset_eff = model.params_additive_offset;
model.params_multiplicative_offset_eff = model.params_multiplicative_offset;
model.params_hyperprior_params_eff = model.params_hyperprior_params;
model.params_family_eff(model.set_default) = [];
model.params_additive_offset_eff(model.set_default) = [];
model.params_multiplicative_offset_eff(model.set_default) = [];
for jj = 1:numel(model.params_hyperprior_params)
    model.params_hyperprior_params_eff{jj}(model.set_default) = [];
end
model.params_theta = rand(1, model.num_eff_params);
model.params_theta(model.params_family_eff == 1) = betainv(model.params_theta(model.params_family_eff == 1), ...
    model.params_hyperprior_params_eff{1}(model.params_family_eff == 1), ...
    model.params_hyperprior_params_eff{2}(model.params_family_eff == 1));
model.params_theta(model.params_family_eff == 2) = gaminv(model.params_theta(model.params_family_eff == 2), ...
    model.params_hyperprior_params_eff{1}(model.params_family_eff == 2), ...
    model.params_hyperprior_params_eff{2}(model.params_family_eff == 2));
model.params_theta(model.params_family_eff == 3) = norminv(model.params_theta(model.params_family_eff == 3), ...
    model.params_hyperprior_params_eff{1}(model.params_family_eff == 3), ...
    sqrt(model.params_hyperprior_params_eff{2}(model.params_family_eff == 3)));
model.params_theta = model.params_multiplicative_offset_eff .* model.params_theta + model.params_additive_offset_eff;
model.params_phi = model.theta_phi(model.params_theta, model);

model.simulate_data = @(model, num_repeats) simulate_data(model, num_repeats);

model.simulate_prob_resp = @(params, model) simulate_prob_resp(params, model);
model.simulate_prob_resp_pred = @(params, model) simulate_prob_resp_pred(params, model);

% model.plot_model_pred=@(model,figid,showlegend) plot_model_pred(model,figid,showlegend);
model.plot_model_pred = @(model, fig_handle, showlegend) plot_model_pred(model, fig_handle, showlegend);

model.get_samples_posterior = @(model, numchains, num_samples, num_samples_save) get_samples_posterior(model, numchains, num_samples, num_samples_save);
model.get_samples_posterior_parallel = @(model, numchains, num_samples, num_samples_save) get_samples_posterior_parallel(model, numchains, num_samples, num_samples_save);

model.log_unnorm_post = @(params, model) log_unnorm_post(params, model);
model.log_likelihood = @(params, model) log_likelihood(params, model);
model.log_prior_theta = @(theta, model) log_prior_theta(theta, model);
model.log_prior_phi = @(phi, model) log_prior_phi(phi, model);

model.get_map = @(model, nstart, gethessian, useparallel) get_map(model, nstart, gethessian, useparallel);
model.get_mle = @(model, nstart, gethessian, useparallel) get_mle(model, nstart, gethessian, useparallel);

model.get_bootstraps_map = @(model, nboot) get_bootstraps_map(model, nboot);
model.get_bootstraps_mle = @(model, nboot) get_bootstraps_mle(model, nboot);

end

function theta = phi_theta(phi, model)
phi=min(model.max_phi,max(-1*model.max_phi,phi));
theta = normcdf(phi, 0, model.normprior_std);
if sum(model.params_family_eff == 1)>0
    theta(:,model.params_family_eff == 1) = betainv(theta(:,model.params_family_eff == 1), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 1),size(phi,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 1),size(phi,1),1));
end
if sum(model.params_family_eff == 2)>0
    theta(:,model.params_family_eff == 2) = gaminv(theta(:,model.params_family_eff == 2), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 2),size(phi,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 2),size(phi,1),1));
end
if sum(model.params_family_eff == 3)>0
    theta(:,model.params_family_eff == 3) = norminv(theta(:,model.params_family_eff == 3), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 3),size(phi,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 3),size(phi,1),1));
end
if sum(model.params_family_eff == 4)>0
    theta(:,model.params_family_eff == 4) = logninv(theta(:,model.params_family_eff == 4), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 4),size(phi,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 4),size(phi,1),1));
end
if sum(model.params_family_eff == 5)>0
    theta(:,model.params_family_eff == 5) = betaprinv(theta(:,model.params_family_eff == 5), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 5),size(phi,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 5),size(phi,1),1));
end



theta = theta.*repmat(model.params_multiplicative_offset_eff,size(theta,1),1) + repmat(model.params_additive_offset_eff,size(theta,1),1);
end

function phi = theta_phi(theta, model)

phi = (theta - repmat(model.params_additive_offset_eff,size(theta,1),1)) ./ repmat(model.params_multiplicative_offset_eff,size(theta,1),1);
if sum(model.params_family_eff == 1)>0
    phi(:,model.params_family_eff == 1) = betacdf(phi(:,model.params_family_eff == 1), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 1),size(phi,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 1),size(phi,1),1));
end
if sum(model.params_family_eff == 2)>0
    phi(:,model.params_family_eff == 2) = gamcdf(phi(:,model.params_family_eff == 2), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 2),size(phi,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 2),size(phi,1),1));
end
if sum(model.params_family_eff == 3)>0
    phi(:,model.params_family_eff == 3) = normcdf(phi(:,model.params_family_eff == 3), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 3),size(phi,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 3),size(phi,1),1));
end
if sum(model.params_family_eff == 4)>0
    phi(:,model.params_family_eff == 4) = logncdf(phi(:,model.params_family_eff == 4), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 4),size(phi,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 4),size(phi,1),1));
end
if sum(model.params_family_eff == 5)>0
    phi(:,model.params_family_eff == 5) = betaprcdf(phi(:,model.params_family_eff == 5), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 5),size(phi,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 5),size(phi,1),1));
end

phi = norminv(phi, 0, model.normprior_std);
phi=min(model.max_phi,max(-1*model.max_phi,phi));
end

function plot_model_pred(model, fig_handle, showlegend)

end


function model = get_samples_posterior(model, numchains, num_samples, num_samples_save)
model = gess_model(model, numchains, num_samples, num_samples_save, 0, 0);
end

function model = get_samples_posterior_parallel(model, numchains, num_samples, num_samples_save)
if isfield(model, 'incomplete_dir_save')
    
    model = gess_parallel_v2(model.incomplete_dir_save,model, numchains, num_samples, num_samples_save, 0, 0);
else
    model = gess_parallel(model, numchains, num_samples, num_samples_save, 0, 0);
end

end

function model1 = get_bootstraps_map(model, nboot)


options = optimoptions(@fminunc, 'Display', 'iter', 'MaxFunctionEvaluations', 1e4,'MaxIterations',1e4);

for i = 1:nboot
    model1{i} = model;
    model1{i}.data{1}.num_ch1 = binornd(model1{i}.data{1}.num_repeats, model1{i}.data{1}.num_ch1 ./ model1{i}.data{1}.num_repeats);
end

parfor i = 1:nboot
    i
    
    param0 = normrnd(0, model.normprior_std, [1, model1{i}.num_eff_params]);
    boots(i, :) = fminunc(@(x) - 1 * (model.log_likelihood(x, model) + model.log_prior_theta(model.phi_theta(x, model), model)), param0, options);
end
model.boots = boots;
end

function model1 = get_bootstraps_mle(model, nboot)


options = optimoptions(@fminunc, 'Display', 'iter', 'MaxFunctionEvaluations', 1e4,'MaxIterations',1e4);

for i = 1:nboot
    model1{i} = model;
    model1{i}.data{1}.num_ch1 = binornd(model1{i}.data{1}.num_repeats, model1{i}.data{1}.num_ch1 ./ model1{i}.data{1}.num_repeats);
end

parfor i = 1:nboot
    i
    
    param0 = normrnd(0, model.normprior_std, [1, model1{i}.num_eff_params]);
    boots(i, :) = fminunc(@(x) - 1 * (model.log_likelihood(x, model)), param0, options);
end
model.boots = boots;
end

function model = get_map(model, nstart, gethessian, useparallel)

options = optimoptions(@fminunc, 'Display', 'iter', 'MaxFunctionEvaluations', 1e4,'MaxIterations',1e4);
if useparallel == 1
    parfor i = 1:nstart
        param0 = normrnd(0, model.normprior_std, [1, model.num_eff_params]);
        if gethessian == 1
            
            [fit(i, :), neg_lups(i), ~, ~, ~, hes{i}] = fminunc(@(x) - 1 * (model.log_likelihood(x, model) + model.log_prior_phi(x, model)), param0, options);
        else
            [fit(i, :), neg_lups(i)] = fminunc(@(x) - 1 * (model.log_likelihood(x, model) + model.log_prior_phi(x, model)), param0, options);
            %             [fit(i,:),neg_lups(i)]=fminunc(@(x) -1*model.log_unnorm_post(x,model),param0,options);
            %         [fit(i,:),neg_lups(i)]=bads(@(x) -1*model.log_unnorm_post(x,model),param0,-100*ones(1,model.num_params),100*ones(1,model.num_params));
        end
    end
else
    
    for i = 1:nstart
        param0 = normrnd(0, model.normprior_std, [1, model.num_eff_params]);
        if gethessian == 1
            [fit(i, :), neg_lups(i), ~, ~, ~, hes{i}] = fminunc(@(x) - 1 * (model.log_likelihood(x, model) + model.log_prior_phi(x, model)), param0, options);
        else
            [fit(i, :), neg_lups(i)] = fminunc(@(x) - 1 * (model.log_likelihood(x, model) + model.log_prior_phi(x, model)), param0, options);
        end
    end
end

model.params_phi_map = fit(find(neg_lups == min(neg_lups), 1), :);
model.map_allparams=fit;
model.map_neg_lups=neg_lups;
if gethessian == 1
    model.params_phi_hessian = hes{find(neg_lups == min(neg_lups), 1)};
    
end
end

function model = get_mle(model, nstart, gethessian, useparallel)

options = optimoptions(@fminunc, 'Display', 'iter', 'MaxFunctionEvaluations', 1e4,'MaxIterations',1e4);
if useparallel == 1
    parfor i = 1:nstart
        param0 = normrnd(0, model.normprior_std, [1, model.num_eff_params]);
        if gethessian == 1
            
            [fit(i, :), neg_lups(i), ~, ~, ~, hes{i}] = fminunc(@(x) - 1 * (model.log_likelihood(x, model)), param0, options);
        else
            [fit(i, :), neg_lups(i)] = fminunc(@(x) - 1 * (model.log_likelihood(x, model)), param0, options);
        end
    end
else
    
    for i = 1:nstart
        param0 = normrnd(0, model.normprior_std, [1, model.num_eff_params]);
        if gethessian == 1
            [fit(i, :), neg_lups(i), ~, ~, ~, hes{i}] = fminunc(@(x) - 1 * (model.log_likelihood(x, model)), param0, options);
        else
            [fit(i, :), neg_lups(i)] = fminunc(@(x) - 1 * (model.log_likelihood(x, model)), param0, options);
        end
    end
end

model.params_phi_mle = fit(find(neg_lups == min(neg_lups), 1), :);
model.mle_allparams=fit;
model.mle_neg_ll=neg_lups;
if gethessian == 1
    model.params_phi_hessian = hes{find(neg_lups == min(neg_lups), 1)};
    
end
end

function lp_theta = log_prior_theta(theta, model)
lp_theta = 0;
theta = (theta - repmat(model.params_additive_offset_eff,size(theta,1),1)) ./ repmat(model.params_multiplicative_offset_eff,size(theta,1),1);

if sum(model.params_family_eff == 1)>0
    lp_theta = lp_theta + sum(log(betapdf(theta(:,model.params_family_eff == 1), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 1),size(theta,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 1),size(theta,1),1))),2);
    
end
if sum(model.params_family_eff == 2)>0
    lp_theta = lp_theta + sum(log(gampdf(theta(:,model.params_family_eff == 2), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 2),size(theta,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 2),size(theta,1),1))),2);
    
end
if sum(model.params_family_eff == 3)>0
    lp_theta = lp_theta + sum(lognormpdf(theta(:,model.params_family_eff == 3), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 3),size(theta,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 3),size(theta,1),1)),2);
    
end

if sum(model.params_family_eff == 4)>0
    lp_theta = lp_theta + sum(log(lognpdf(theta(:,model.params_family_eff == 4), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 4),size(theta,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 4),size(theta,1),1))),2);
    
    
end

if sum(model.params_family_eff == 5)>0
    lp_theta = lp_theta + sum(logbetaprpdf(theta(:,model.params_family_eff == 5), ...
        repmat(model.params_hyperprior_params_eff{1}(model.params_family_eff == 5),size(theta,1),1), ...
        repmat(model.params_hyperprior_params_eff{2}(model.params_family_eff == 5),size(theta,1),1)),2);
    
end


if isinf(lp_theta) || isnan(lp_theta)
    warning('inf or nan log prior');
    lp_theta=-1e4;
    %keyboard
end
lp_theta= lp_theta - sum(log(repmat(model.params_multiplicative_offset_eff,size(theta,1),1)),2);

end

function lp_phi = log_prior_phi(phi, model)


lp_phi=sum(lognormpdf(phi,0,model.normprior_std),2);


end

function ll = log_likelihood(params_phi, model)

p_ch=simulate_prob_resp(model.phi_theta(params_phi,model),model);
ll=0;
for jj=1:numel(model.design_matrix)
    if any(p_ch{jj}<=0) || any(p_ch{jj}>=1)
        warning('0 or 1 probability of choice');
        %         keyboard
        p_ch{jj}=min(1-eps,max(eps,p_ch{jj}));
    end
    
    ll=ll+sum((model.data{jj}.num_ch1.*log(p_ch{jj}))+((model.data{jj}.num_repeats-model.data{jj}.num_ch1).*log(1-p_ch{jj})),2);
    
    
end

end
function lup = log_unnorm_post(params_phi, model)
lup = model.log_likelihood(params_phi, model) + model.log_prior_phi(params_phi,model);
end

function model = simulate_data(model, num_repeats)
p_ch1=simulate_prob_resp(model.phi_theta(model.params_phi,model),model);
for jj=1:numel(model.design_matrix)
    model.data{jj}.num_repeats=num_repeats*ones(1,model.design_matrix{jj}.npts);
    model.data{jj}.num_ch1=binornd(num_repeats,p_ch1{jj});
end
end

function theta = get_full_params(theta0, model)

theta = zeros(1, model.num_params);

theta(model.set_default) = model.default_values;
theta(setdiff([1:model.num_params], model.set_default)) = theta0;

end

function p_ch1 = simulate_prob_resp(params, model)

params=get_full_params(params,model);





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



end



function pr=pr_rl_il(pi,sigmap_2,sigmaf_2,Il)

rls=(dec2bin(0:2^3-1)-'0')+1;
rls(:,1)=rls(:,1)+1;
rls(:,2)=rls(:,2)+3;
rls(:,3)=rls(:,3)+5;
for i=1:size(rls,1)
    
end

end



function pr=pr_rt_il_rl(pi,sigmap_2,sigmaf_2,Il,rl)
Il(rl==1)=Il(rl==1)/sigmaf_2;
Il(rl==0)=Il(rl==0)/sigmap_2;
pr=sigmoid(exp(sum(log(pi/(1-pi))+log(sigmoid(2*Il-log(pi/(1-pi))))-log(sigmoid(2*Il+log(pi/(1-pi)))))));
end

