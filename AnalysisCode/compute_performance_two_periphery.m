function [mn_sig,mn_sig_norm,mn_sig_ideal_foveal_only,mn_sig_ideal_foveal_only_norm,...
    mn_sig_ideal_all,mn_sig_ideal_all_norm,mn_perf,mn_perf_ideal_foveal_only,mn_perf_ideal_all,...
    tr_seg,err_perf,tr_seg_ideal_foveal_only,err_perf_ideal_foveal_only,tr_seg_ideal_all,err_perf_ideal_all]...
    = compute_performance_two_periphery(temporal_kernel,signal_all,signal_chosen_raw,...
    accuracy,accuracy_ideal_all,accuracy_ideal_foveal_only,all_frames,segments)

% segments = 3;
mn_sig = [];
mn_sig_ideal_foveal_only = [];
mn_sig_ideal_all = [];
mn_perf = [];
mn_perf_ideal_foveal_only = [];
mn_perf_ideal_all = [];
tr_seg = [];
err_perf = [];
tr_seg_ideal_foveal_only = [];
err_perf_ideal_foveal_only = [];
tr_seg_ideal_all = [];
err_perf_ideal_all = [];

weighted_signal_all = abs(sum(temporal_kernel(1:all_frames) .* signal_all,2));
[sorted_weighted_signal_all,order] = sort(weighted_signal_all);
sorted_accuracy = accuracy(order);

bin_edges = linspace(sorted_weighted_signal_all(1),sorted_weighted_signal_all(end),segments+1);
for j=1:segments
    ind_chosen = find(sorted_weighted_signal_all>=bin_edges(j) & sorted_weighted_signal_all<=bin_edges(j+1));
    mn_sig(j) = (bin_edges(j) + bin_edges(j+1))/2;
    mn_perf(j) = sum(sorted_accuracy(ind_chosen))/length(ind_chosen);
    err_perf(j) = sqrt(mn_perf(j)*(1-mn_perf(j))/length(ind_chosen));
    tr_seg(j) = length(ind_chosen);
end
mn_sig_norm = mn_sig/mn_sig(end);


signal_ideal_all = abs(sum(signal_all,2));
[sorted_signal_ideal_all,order1] = sort(signal_ideal_all);
sorted_accuracy_ideal_all = accuracy_ideal_all(order1);
bin_edges_ideal_all = linspace(sorted_signal_ideal_all(1),sorted_signal_ideal_all(end),segments+1);

for j=1:segments
    ind_chosen2 = find(sorted_signal_ideal_all>=bin_edges_ideal_all(j) & sorted_signal_ideal_all<=bin_edges_ideal_all(j+1));
    mn_sig_ideal_all(j) = (bin_edges_ideal_all(j) + bin_edges_ideal_all(j+1))/2;
    mn_perf_ideal_all(j) = sum(sorted_accuracy_ideal_all(ind_chosen2))/length(ind_chosen2);
    err_perf_ideal_all(j) = sqrt(mn_perf_ideal_all(j)*(1-mn_perf_ideal_all(j))/length(ind_chosen2));
    tr_seg_ideal_all(j) = length(ind_chosen2);
end
mn_sig_ideal_all_norm = mn_sig_ideal_all/mn_sig_ideal_all(end);


signal_ideal_foveal_only = abs(sum(signal_chosen_raw,2));
[sorted_signal_ideal_foveal_only,order2] = sort(signal_ideal_foveal_only);
sorted_accuracy_ideal_foveal_only = accuracy_ideal_foveal_only(order2);
bin_edges_ideal_foveal_only = linspace(sorted_signal_ideal_foveal_only(1),sorted_signal_ideal_foveal_only(end),segments+1);

for j=1:segments
    ind_chosen1 = find(sorted_signal_ideal_foveal_only>=bin_edges_ideal_foveal_only(j) & sorted_signal_ideal_foveal_only<=bin_edges_ideal_foveal_only(j+1));
    mn_sig_ideal_foveal_only(j) = (bin_edges_ideal_foveal_only(j) + bin_edges_ideal_foveal_only(j+1))/2;
    mn_perf_ideal_foveal_only(j) = sum(sorted_accuracy_ideal_foveal_only(ind_chosen1))/length(ind_chosen1);
    err_perf_ideal_foveal_only(j) = sqrt(mn_perf_ideal_foveal_only(j)*(1-mn_perf_ideal_foveal_only(j))/length(ind_chosen1));
    tr_seg_ideal_foveal_only(j) = length(ind_chosen1);
end
mn_sig_ideal_foveal_only_norm = mn_sig_ideal_foveal_only/mn_sig_ideal_foveal_only(end);
