%Eq 6

function [score, comp1, comp2] = compute_score_vm_v2(table, l_dim, num_samples, is_sampling)
% rng(seed)
posterior_c_t = marginalize_table(table, 2:3);
if is_sampling==1
    num_samples_neg = binornd(num_samples, posterior_c_t(1));
    p_neg = num_samples_neg / num_samples;
    p_pos = 1 - p_neg;
else
    p_neg = posterior_c_t(1);
    p_pos = posterior_c_t(2);
end

conditional1 = normalize(table(1, :, :));
conditional2 = normalize(table(2, :, :));
% table_sampled = cat(1, p_neg * marginal1, p_pos * marginal2);
table_sampled=p_neg*conditional1 + p_pos*conditional2;
% table_sampled=posterior_c_t(1)*conditional1 + posterior_c_t(2)*conditional2;

comp1 = entropy_compute(table_sampled, 1);


if l_dim==2
    posterior_ldim = marginalize_table(table, [1,3]);
    
    comp2=0;
    for j=1:length(posterior_ldim)
        
%         cond_posterior_c_t = marginalize_table(squeeze(normalize(table(:, j, :))),2);
%         if is_sampling==1
%             num_samples_neg = binornd(num_samples, cond_posterior_c_t(1));
%             p_neg = num_samples_neg / num_samples;
%             p_pos = 1 - p_neg;
%         else
%             p_neg = cond_posterior_c_t(1);
%             p_pos = cond_posterior_c_t(2);
%         end
%         marginal1 = normalize(table(1, j, :));
%         marginal2 = normalize(table(2, j, :));
%         
%         table_sampled=p_neg*marginal1 + p_pos*marginal2;
        
        
%         comp2 = comp2 + posterior_ldim(j)*entropy_compute(table_sampled, 1);
        comp2 = comp2 + posterior_ldim(j)*entropy_compute(normalize(table_sampled(:,j,:)), 1);
        
        
    end
else
    posterior_ldim = marginalize_table(table, [1,2]);
    
    comp2=0;
    for j=1:length(posterior_ldim)
        
%         cond_posterior_c_t = marginalize_table(squeeze(normalize(table(:, :, j))),2);
%         if is_sampling==1
%             num_samples_neg = binornd(num_samples, cond_posterior_c_t(1));
%             p_neg = num_samples_neg / num_samples;
%             p_pos = 1 - p_neg;
%         else
%             p_neg = cond_posterior_c_t(1);
%             p_pos = cond_posterior_c_t(2);
%         end
%         marginal1 = normalize(table(1, :, j));
%         marginal2 = normalize(table(2, :, j));
%         
%         table_sampled=p_neg*marginal1 + p_pos*marginal2;
%         
%         
%         comp2 = comp2 + posterior_ldim(j)*entropy_compute(table_sampled, 1);
        comp2 = comp2 + posterior_ldim(j)*entropy_compute(normalize(table_sampled(:,:,j)), 1);
        
        
    end
    
end


score = comp1 - comp2;
end

function tbl = normalize(tbl)
tbl = tbl / sum(tbl(:));
end