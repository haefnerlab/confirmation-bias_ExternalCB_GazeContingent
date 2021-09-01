function table_ct_cl1_cl2_given_d = compute_table(I_peri_l1,I_peri_l2,sigma_peri,prob_ct,p_match,log_odds_accumulated)
num_dim = 3;
dims=[2,2,2];
log_table_ct_cl1_cl2_given_d = zeros(dims);
index_val_array = dec2bin(2^num_dim-1:-1:0)-'0';
val_array = index_val_array;
val_array = double(val_array);
val_array(val_array==0) = -1;
idx_ct = 1;
idx_cl1 = 2;
idx_cl2 = 3;
loglik_ct=log([sigmoid(log_odds_accumulated),1-sigmoid(log_odds_accumulated)]);

for i=1:size(index_val_array,1)
    log_temp = loglik_ct(index_val_array(i,idx_ct)+1) ...
        + lognormpdf(I_peri_l1,val_array(i,idx_cl1),sigma_peri) ...
        + lognormpdf(I_peri_l2,val_array(i,idx_cl2),sigma_peri) ...
        + log(prob_ct(index_val_array(i,idx_ct)+1)) ...
        + log(p_match^(val_array(i,idx_cl1)==val_array(i,idx_ct)) * (1-p_match)^((val_array(i,idx_cl1)~=val_array(i,idx_ct)))) ...
        + log(p_match^(val_array(i,idx_cl2)==val_array(i,idx_ct)) * (1-p_match)^((val_array(i,idx_cl2)~=val_array(i,idx_ct))));
   
    log_table_ct_cl1_cl2_given_d(...
       index_val_array(i,idx_ct)+1,...
       index_val_array(i,idx_cl1)+1,...
       index_val_array(i,idx_cl2)+1) = log_temp;
end
table_ct_cl1_cl2_given_d = exp(log_table_ct_cl1_cl2_given_d-logsumexp2(log_table_ct_cl1_cl2_given_d(:)'));
end