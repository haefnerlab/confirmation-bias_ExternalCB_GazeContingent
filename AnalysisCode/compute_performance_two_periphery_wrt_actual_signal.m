function [mn_sig_ideal_all,mn_perf,mn_perf_ideal_foveal_only,mn_perf_ideal_all,...
    tr_seg,err_perf,tr_seg_ideal_foveal_only,err_perf_ideal_foveal_only,tr_seg_ideal_all,err_perf_ideal_all]...
    = compute_performance_two_periphery_wrt_actual_signal(signal_all,...
    accuracy,accuracy_ideal_all,accuracy_ideal_foveal_only,segments)

signal_ideal_all = abs(sum(signal_all,2));
[sorted_signal_ideal_all,order] = sort(signal_ideal_all);
sorted_accuracy_ideal_all = accuracy_ideal_all(order);
sorted_accuracy = accuracy(order);
sorted_accuracy_ideal_foveal_only = accuracy_ideal_foveal_only(order);

bin_edges_ideal_all = linspace(sorted_signal_ideal_all(1),sorted_signal_ideal_all(end),segments+1);

for j=1:segments
    ind_chosen = find(sorted_signal_ideal_all>=bin_edges_ideal_all(j) & sorted_signal_ideal_all<bin_edges_ideal_all(j+1));
    mn_sig_ideal_all(j) = (bin_edges_ideal_all(j) + bin_edges_ideal_all(j+1))/2;
    
    mn_perf_ideal_all(j) = sum(sorted_accuracy_ideal_all(ind_chosen))/length(ind_chosen);
    err_perf_ideal_all(j) = sqrt(mn_perf_ideal_all(j)*(1-mn_perf_ideal_all(j))/length(ind_chosen));
    tr_seg_ideal_all(j) = length(ind_chosen);
    
    mn_perf(j) = sum(sorted_accuracy(ind_chosen))/length(ind_chosen);
    err_perf(j) = sqrt(mn_perf(j)*(1-mn_perf(j))/length(ind_chosen));
    tr_seg(j) = length(ind_chosen);
    
    mn_perf_ideal_foveal_only(j) = sum(sorted_accuracy_ideal_foveal_only(ind_chosen))/length(ind_chosen);
    err_perf_ideal_foveal_only(j) = sqrt(mn_perf_ideal_foveal_only(j)*(1-mn_perf_ideal_foveal_only(j))/length(ind_chosen));
    tr_seg_ideal_foveal_only(j) = length(ind_chosen);
end
