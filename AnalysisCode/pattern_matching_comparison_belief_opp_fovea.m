function [one_previous, one_previous_belief]...
    = pattern_matching_comparison_belief_opp_fovea(data, temporal_kernel, bias)

frames = data.number_of_images;
num_peri = 2;
signal_chosen_raw_actual = [];
signal_not_chosen_raw_actual = [];
signal_all_actual = [];
choice_raw = [];
% we store the signals w.r.t the number of elements on the periphery
for k = 1:data.current_trial-1
    if data.num_peri(k)==2
        [chosen_ideal_frame_signals,notchosen_ideal_frame_signals,~] = regenerate_signals(data,k);
        signal_chosen_raw_actual = [signal_chosen_raw_actual; chosen_ideal_frame_signals'];
        signal_not_chosen_raw_actual = [signal_not_chosen_raw_actual; notchosen_ideal_frame_signals'];
        signal_all_actual = [signal_all_actual; chosen_ideal_frame_signals' notchosen_ideal_frame_signals'];
        choice_raw = [choice_raw data.choice(k)];
    end
end

trials = size(choice_raw, 2);
all_frames = 1 + frames * num_peri;

signal_chosen_raw = sign(signal_chosen_raw_actual);
signal_not_chosen_raw = sign(signal_not_chosen_raw_actual);
signal_all = sign(signal_all_actual);
weighted_signal_all = temporal_kernel .* signal_all;
weighted_signal_chosen = weighted_signal_all (:,1:frames);
weighted_signal_notchosen = weighted_signal_all(:,frames+1:end-num_peri);

evi1 = bias + weighted_signal_chosen(:,1);
evi2 = bias + sum(weighted_signal_chosen(:,1:2),2) + weighted_signal_notchosen(:,1);
evi3 = bias + sum(weighted_signal_chosen(:,1:3),2) + sum(weighted_signal_notchosen(:,1:2),2);

% one previous case
opp_tr = squeeze(evi1 .* signal_all(:,1));
tr = signal_chosen_raw(:,2) .* signal_not_chosen_raw(:,1);

one_previous1 = squeeze(signal_all((tr==-1 & sign(opp_tr)==-1),1) .* signal_all((tr==-1 & sign(opp_tr)==-1),2));
one_previous1(one_previous1==-1) = 0;
one_previous(1,1) = mean(one_previous1);
one_previous(1,2) = sqrt(var(one_previous1)/length(one_previous1));

bf_match1 = sign(squeeze(signal_all((tr==-1 & sign(opp_tr)==-1),2)))==sign(evi1((tr==-1 & sign(opp_tr)==-1)));
one_previous_belief(1,1) = mean(bf_match1);
one_previous_belief(1,2) = sqrt(var(bf_match1)/length(bf_match1));

tr = [];
opp_tr = squeeze(evi2 .* signal_all(:,2));
tr = signal_chosen_raw(:,3) .* signal_not_chosen_raw(:,2);

one_previous2 = squeeze(signal_all((tr==-1 & sign(opp_tr)==-1),2) .* signal_all((tr==-1 & sign(opp_tr)==-1),3));
one_previous2(one_previous2==-1) = 0;
one_previous(2,1) = mean(one_previous2);
one_previous(2,2) = sqrt(var(one_previous2)/length(one_previous2));

bf_match2 = sign(squeeze(signal_all((tr==-1 & sign(opp_tr)==-1 ),3)))==sign(evi2((tr==-1 & sign(opp_tr)==-1)));
one_previous_belief(2,1) = mean(bf_match2);
one_previous_belief(2,2) = sqrt(var(bf_match2)/length(bf_match2));

tr = [];
opp_tr = squeeze(evi3 .* signal_all(:,3));
tr = signal_chosen_raw(:,4) .* signal_not_chosen_raw(:,3);

one_previous3 = squeeze(signal_all((tr==-1 & sign(opp_tr)==-1),3) .* signal_all((tr==-1 & sign(opp_tr)==-1),4));
one_previous3(one_previous3==-1) = 0;
one_previous(3,1) = mean(one_previous3);
one_previous(3,2) = sqrt(var(one_previous3)/length(one_previous3));

bf_match3 = sign(squeeze(signal_all((tr==-1 & sign(opp_tr)==-1),4)))==sign(evi3((tr==-1 & sign(opp_tr)==-1)));
one_previous_belief(3,1) = mean(bf_match3);
one_previous_belief(3,2) = sqrt(var(bf_match3)/length(bf_match3));

one_previous_all = [one_previous1(:); one_previous2(:); one_previous3(:)];
one_previous(4,1) = mean(one_previous_all);
one_previous(4,2) = sqrt(var(one_previous_all)/length(one_previous_all));
bf_match_all = [bf_match1(:); bf_match2(:); bf_match3(:)];
one_previous_belief(4,1) = mean(bf_match_all);
one_previous_belief(4,2) = sqrt(var(bf_match_all)/length(bf_match_all));


end