function [score, comp1, comp2] = compute_score(table, l_dim, num_samples_total_entropy,num_samples_noise_entropy, is_sampling,seed)
rng(seed)
posterior_c_t = marginalize_table(table, 2:3);
if is_sampling(1)==1
    
    
    if rand()<(num_samples_total_entropy-floor(num_samples_total_entropy))
        num_samples_neg_total_entropy = binornd(ceil(num_samples_total_entropy), posterior_c_t(1));
    else
        num_samples_neg_total_entropy = binornd(floor(num_samples_total_entropy), posterior_c_t(1));
    end
    p_neg_total_entropy = num_samples_neg_total_entropy / num_samples_total_entropy;
    p_pos_total_entropy = 1 - p_neg_total_entropy;
    
    
    
    
else
    p_neg_total_entropy = posterior_c_t(1);
    p_pos_total_entropy = posterior_c_t(2);
    
end


if is_sampling(2)==1
    
    
    
    
    if rand()<(num_samples_noise_entropy-floor(num_samples_noise_entropy))
        num_samples_neg_noise_entropy = binornd(ceil(num_samples_noise_entropy), posterior_c_t(1));
    else
        num_samples_neg_noise_entropy = binornd(floor(num_samples_noise_entropy), posterior_c_t(1));
    end
    p_neg_noise_entropy = num_samples_neg_noise_entropy / num_samples_noise_entropy;
    p_pos_noise_entropy = 1 - p_neg_noise_entropy;
    
    
else
    
    p_neg_noise_entropy = posterior_c_t(1);
    p_pos_noise_entropy = posterior_c_t(2);
end


conditional1 = normalize(table(1, :, :));
conditional2 = normalize(table(2, :, :));
% table_sampled = cat(1, p_neg * marginal1, p_pos * marginal2);
table_sampled=p_neg_total_entropy*conditional1 + p_pos_total_entropy*conditional2;
% table_sampled=posterior_c_t(1)*marginal1 + posterior_c_t(2)*marginal2;

comp1 = entropy_compute(table_sampled, l_dim);
% comp1 = entropy_compute(table, l_dim); %%old version

% DEBUG - resample p_neg and p_pos again here
% if is_sampling==1
%     num_samples_neg = binornd(num_samples, posterior_c_t(1));
%     p_neg = num_samples_neg / num_samples;
%     p_pos = 1 - p_neg;
% else
%     p_neg = posterior_c_t(1);
%     p_pos = posterior_c_t(2);
% end

comp2 = p_neg_noise_entropy * entropy_compute(conditional1, l_dim) + ...
    p_pos_noise_entropy * entropy_compute(conditional2, l_dim);

score = comp1 - comp2;
end

function tbl = normalize(tbl)
tbl = tbl / sum(tbl(:));
end