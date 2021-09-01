function [score, comp1, comp2] = compute_score_im(table, l_dim, num_samples, is_sampling)

posterior_c_t = marginalize_table(table, 2:4);
if is_sampling==1
    num_samples_neg = binornd(num_samples, posterior_c_t(1));
    p_neg = num_samples_neg / num_samples;
    p_pos = 1 - p_neg;
else
    p_neg = posterior_c_t(1);
    p_pos = posterior_c_t(2);
end

marginal1 = normalize(table(1, :, :, :));
marginal2 = normalize(table(2, :, :, :));
% table_sampled = cat(1, p_neg * marginal1, p_pos * marginal2);
table_sampled=p_neg*marginal1 + p_pos*marginal2;
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

comp2 = p_neg * entropy_compute(marginal1, l_dim) + ...
    p_pos * entropy_compute(marginal2, l_dim);

score = comp1 - comp2;
end

function tbl = normalize(tbl)
tbl = tbl / sum(tbl(:));
end