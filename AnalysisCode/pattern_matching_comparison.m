function [one_previous, one_previous_belief,...
    two_previous, two_previous_belief,...
    three_previous, three_previous_belief]...
    = pattern_matching_comparison(data, temporal_kernel, bias)

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
tr = signal_chosen_raw(:,2) .* signal_not_chosen_raw(:,1);
one_previous1 = squeeze(signal_all((tr==-1),1) .* signal_all((tr==-1),2));
one_previous1(one_previous1==-1) = 0;
one_previous(1,1) = mean(one_previous1);
one_previous(1,2) = sqrt(var(one_previous1)/length(one_previous1));

bf_match1 = sign(squeeze(signal_all((tr==-1),2)))==sign(evi1((tr==-1)));
one_previous_belief(1,1) = mean(bf_match1);
one_previous_belief(1,2) = sqrt(var(bf_match1)/length(bf_match1));

tr = [];
tr = signal_chosen_raw(:,3) .* signal_not_chosen_raw(:,2);
one_previous2 = squeeze(signal_all((tr==-1),2) .* signal_all((tr==-1),3));
one_previous2(one_previous2==-1) = 0;
one_previous(2,1) = mean(one_previous2);
one_previous(2,2) = sqrt(var(one_previous2)/length(one_previous2));

bf_match2 = sign(squeeze(signal_all((tr==-1),3)))==sign(evi2((tr==-1)));
one_previous_belief(2,1) = mean(bf_match2);
one_previous_belief(2,2) = sqrt(var(bf_match2)/length(bf_match2));

tr = [];
tr = signal_chosen_raw(:,4) .* signal_not_chosen_raw(:,3);
one_previous3 = squeeze(signal_all((tr==-1),3) .* signal_all((tr==-1),4));
one_previous3(one_previous3==-1) = 0;
one_previous(3,1) = mean(one_previous3);
one_previous(3,2) = sqrt(var(one_previous3)/length(one_previous3));

bf_match3 = sign(squeeze(signal_all((tr==-1),4)))==sign(evi3((tr==-1)));
one_previous_belief(3,1) = mean(bf_match3);
one_previous_belief(3,2) = sqrt(var(bf_match3)/length(bf_match3));

one_previous_all = [one_previous1(:); one_previous2(:); one_previous3(:)];
one_previous(4,1) = mean(one_previous_all);
one_previous(4,2) = sqrt(var(one_previous_all)/length(one_previous_all));
bf_match_all = [bf_match1(:); bf_match2(:); bf_match3(:)];
one_previous_belief(4,1) = mean(bf_match_all);
one_previous_belief(4,2) = sqrt(var(bf_match_all)/length(bf_match_all));


% two previous case
tr = [];
tr = signal_chosen_raw(:,3) .* signal_not_chosen_raw(:,2); 
temp = squeeze(signal_all(:,1) .* signal_all(:,2));
two_previous2_same = squeeze(signal_all((temp==1 & tr==-1)',2) .* signal_all((temp==1 & tr==-1)',3));
two_previous2_same(two_previous2_same==-1) = 0;
two_previous(1,1,1) = mean(two_previous2_same);
two_previous(1,1,2) = sqrt(var(two_previous2_same)/length(two_previous2_same));

bf_match4 = sign(squeeze(signal_all((temp==1 & tr==-1),3)))==sign(evi2((temp==1 & tr==-1)));
two_previous_belief(1,1,1) = mean(bf_match4);
two_previous_belief(1,1,2) = sqrt(var(bf_match4)/length(bf_match4));

two_previous2_diff = squeeze(signal_all((temp==-1 & tr==-1)',2) .* signal_all((temp==-1 & tr==-1)',3));
two_previous2_diff(two_previous2_diff==-1) = 0;
two_previous(1,2,1) = mean(two_previous2_diff);
two_previous(1,2,2) = sqrt(var(two_previous2_diff)/length(two_previous2_diff));

bf_match5 = sign(squeeze(signal_all((temp==-1 & tr==-1),3)))==sign(evi2((temp==-1 & tr==-1)));
two_previous_belief(1,2,1) = mean(bf_match5);
two_previous_belief(1,2,2) = sqrt(var(bf_match5)/length(bf_match5));

two_previous2_all = [two_previous2_same(:); two_previous2_diff(:)];
two_previous(1,3,1) = mean(two_previous2_all);
two_previous(1,3,2) = sqrt(var(two_previous2_all)/length(two_previous2_all));

bf_match6 = [bf_match4(:); bf_match5(:)];
two_previous_belief(1,3,1) = mean(bf_match6);
two_previous_belief(1,3,2) = sqrt(var(bf_match6)/length(bf_match6));


temp = [];
tr = [];
tr = signal_chosen_raw(:,4) .* signal_not_chosen_raw(:,3);
temp = squeeze(signal_all(:,2) .* signal_all(:,3));
two_previous3_same = squeeze(signal_all((temp==1 & tr==-1)',3) .* signal_all((temp==1 & tr==-1)',4));
two_previous3_same(two_previous3_same==-1) = 0;
two_previous(2,1,1) = mean(two_previous3_same);
two_previous(2,1,2) = sqrt(var(two_previous3_same)/length(two_previous3_same));

bf_match7 = sign(squeeze(signal_all((temp==1 & tr==-1),4)))==sign(evi3((temp==1 & tr==-1)));
two_previous_belief(2,1,1) = mean(bf_match7);
two_previous_belief(2,1,2) = sqrt(var(bf_match7)/length(bf_match7));

two_previous3_diff = squeeze(signal_all((temp==-1 & tr==-1)',3) .* signal_all((temp==-1 & tr==-1)',4));
two_previous3_diff(two_previous3_diff==-1) = 0;
two_previous(2,2,1) = mean(two_previous3_diff);
two_previous(2,2,2) = sqrt(var(two_previous3_diff)/length(two_previous3_diff));

bf_match8 = sign(squeeze(signal_all((temp==-1 & tr==-1),4)))==sign(evi3((temp==-1 & tr==-1)));
two_previous_belief(2,2,1) = mean(bf_match8);
two_previous_belief(2,2,2) = sqrt(var(bf_match8)/length(bf_match8));

two_previous3_all = [two_previous3_same(:); two_previous3_diff(:)];
two_previous(2,3,1) = mean(two_previous3_all);
two_previous(2,3,2) = sqrt(var(two_previous3_all)/length(two_previous3_all));

bf_match9 = [bf_match7(:); bf_match8(:)];
two_previous_belief(2,3,1) = mean(bf_match9);
two_previous_belief(2,3,2) = sqrt(var(bf_match9)/length(bf_match9));

two_previous_all_same = [two_previous2_same(:); two_previous3_same(:)];
two_previous(3,1,1) = mean(two_previous_all_same);
two_previous(3,1,2) = sqrt(var(two_previous_all_same)/length(two_previous_all_same));

bf_match10 = [bf_match4(:); bf_match7(:)];
two_previous_belief(3,1,1) = mean(bf_match10);
two_previous_belief(3,1,2) = sqrt(var(bf_match10)/length(bf_match10));

two_previous_all_diff = [two_previous2_diff(:); two_previous3_diff(:)];
two_previous(3,2,1) = mean(two_previous_all_diff);
two_previous(3,2,2) = sqrt(var(two_previous_all_diff)/length(two_previous_all_diff));

bf_match11 = [bf_match5(:); bf_match8(:)];
two_previous_belief(3,2,1) = mean(bf_match11);
two_previous_belief(3,2,2) = sqrt(var(bf_match11)/length(bf_match11));

two_previous_all = [two_previous2_all(:); two_previous3_all(:)];
two_previous(3,3,1) = mean(two_previous_all);
two_previous(3,3,2) = sqrt(var(two_previous_all)/length(two_previous_all));

bf_match12 = [bf_match6(:); bf_match9(:)];
two_previous_belief(3,3,1) = mean(bf_match12);
two_previous_belief(3,3,2) = sqrt(var(bf_match12)/length(bf_match12));


% three previous case
temp = [];
tr = [];
tr = signal_chosen_raw(:,4) .* signal_not_chosen_raw(:,3);
temp = abs(squeeze(signal_all(:,1) + signal_all(:,2) + signal_all(:,3)));
three_previous3_same = squeeze(signal_all((temp==3 & tr==-1)',3) .* signal_all((temp==3 & tr==-1)',4));
three_previous3_same(three_previous3_same==-1) = 0;
three_previous(1,1) = mean(three_previous3_same);
three_previous(1,2) = sqrt(var(three_previous3_same)/length(three_previous3_same));

bf_match13 = sign(squeeze(signal_all((temp==3 & tr==-1),4)))==sign(evi3((temp==3 & tr==-1)));
three_previous_belief(1,1) = mean(bf_match13);
three_previous_belief(1,2) = sqrt(var(bf_match13)/length(bf_match13));

three_previous3_diff = squeeze(signal_all((temp==1 & tr==-1)',3) .* signal_all((temp==1 & tr==-1)',4));
three_previous3_diff(three_previous3_diff==-1) = 0;
three_previous(2,1) = mean(three_previous3_diff);
three_previous(2,2) = sqrt(var(three_previous3_diff)/length(three_previous3_diff));

bf_match14 = sign(squeeze(signal_all((temp==1 & tr==-1),4)))==sign(evi3((temp==1 & tr==-1)));
three_previous_belief(2,1) = mean(bf_match14);
three_previous_belief(2,2) = sqrt(var(bf_match14)/length(bf_match14));

three_previous3_all = [three_previous3_same(:); three_previous3_diff(:)];
three_previous(3,1) = mean(three_previous3_all);
three_previous(3,2) = sqrt(var(three_previous3_all)/length(three_previous3_all));

bf_match_15 = [bf_match13(:); bf_match14(:)];
three_previous_belief(3,1) = mean(bf_match_15);
three_previous_belief(3,2) = sqrt(var(bf_match_15)/length(bf_match_15));

end