function [actual_data, predicted_data, mn_sig] = check_weighted_prediction_vs_true_prediction(choice, signal_all, weighted_signal_all, bias, lapse, bins)
trials = size(signal_all);
for i=1:trials
    summed_sig(i) = sum(signal_all(i,:));
    pred_prob(i) = lapse * 0.5 + (1 - lapse) * sigmoid(bias + sum(weighted_signal_all(i,:)));
end
[sorted_summed_sig, order] = sort(summed_sig);
sorted_choice = choice(order);
sorted_pred_prob = pred_prob(order);
bin_edges = linspace(sorted_summed_sig(1),sorted_summed_sig(end),bins+1);
for j=1:bins
    ind_chosen = find(sorted_summed_sig>=bin_edges(j) & sorted_summed_sig<=bin_edges(j+1));
    mn_sig(j) = (bin_edges(j) + bin_edges(j+1))/2;
    actual_data(j,1) = sum(sorted_choice(ind_chosen))/length(ind_chosen);
    actual_data(j,2) = sqrt(actual_data(j,1)*(1-actual_data(j,1))/length(ind_chosen));
    predicted_data(j,1) = sum(sorted_pred_prob(ind_chosen))/length(ind_chosen);
    predicted_data(j,2) = sqrt(var(sorted_choice(ind_chosen))/length(ind_chosen));
end
end