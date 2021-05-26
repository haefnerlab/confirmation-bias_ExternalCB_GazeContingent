function [one_previous, two_previous, three_previous] = pattern_matching(subjectID, expt_type)

% initialize variables
[ParentFolderPath] = fileparts(pwd);
datadir = fullfile(ParentFolderPath, '/RawData');

signal_chosen_raw_actual = [];
signal_not_chosen_raw_actual = [];
signal_all_actual = [];
choice_raw = [];

% load data
[data,~] = LoadAllSubjectData(subjectID,expt_type,datadir);
frames = data.number_of_images;
num_peri = 2;
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

signal_chosen_raw = sign(signal_chosen_raw_actual);
signal_not_chosen_raw = sign(signal_not_chosen_raw_actual);
signal_all = sign(signal_all_actual);

trials = size(choice_raw, 2);
all_frames = 1 + frames * num_peri;

% one previous case
tr = signal_chosen_raw(:,2) .* signal_not_chosen_raw(:,1);
one_previous1 = squeeze(signal_all((tr==-1),1) .* signal_all((tr==-1),2));
one_previous1(one_previous1==-1) = 0;
one_previous(1,1) = mean(one_previous1);
one_previous(1,2) = sqrt(var(one_previous1)/length(one_previous1));
tr = [];
tr = signal_chosen_raw(:,3) .* signal_not_chosen_raw(:,2);
one_previous2 = squeeze(signal_all((tr==-1),2) .* signal_all((tr==-1),3));
one_previous2(one_previous2==-1) = 0;
one_previous(2,1) = mean(one_previous2);
one_previous(2,2) = sqrt(var(one_previous2)/length(one_previous2));
tr = [];
tr = signal_chosen_raw(:,4) .* signal_not_chosen_raw(:,3);
one_previous3 = squeeze(signal_all((tr==-1),3) .* signal_all((tr==-1),4));
one_previous3(one_previous3==-1) = 0;
one_previous(3,1) = mean(one_previous3);
one_previous(3,2) = sqrt(var(one_previous3)/length(one_previous3));
one_previous_all = [one_previous1(:); one_previous2(:); one_previous3(:)];
one_previous(4,1) = mean(one_previous_all);
one_previous(4,2) = sqrt(var(one_previous_all)/length(one_previous_all));

% two previous case
tr = [];
tr = signal_chosen_raw(:,3) .* signal_not_chosen_raw(:,2); 
temp = squeeze(signal_all(:,1) .* signal_all(:,2));
two_previous2_same = squeeze(signal_all((temp==1 & tr==-1)',2) .* signal_all((temp==1 & tr==-1)',3));
two_previous2_same(two_previous2_same==-1) = 0;
two_previous(1,1,1) = mean(two_previous2_same);
two_previous(1,1,2) = sqrt(var(two_previous2_same)/length(two_previous2_same));

two_previous2_diff = squeeze(signal_all((temp==-1 & tr==-1)',2) .* signal_all((temp==-1 & tr==-1)',3));
two_previous2_diff(two_previous2_diff==-1) = 0;
two_previous(1,2,1) = mean(two_previous2_diff);
two_previous(1,2,2) = sqrt(var(two_previous2_diff)/length(two_previous2_diff));

two_previous2_all = [two_previous2_same(:); two_previous2_diff(:)];
two_previous(1,3,1) = mean(two_previous2_all);
two_previous(1,3,2) = sqrt(var(two_previous2_all)/length(two_previous2_all));


temp = [];
tr = [];
tr = signal_chosen_raw(:,4) .* signal_not_chosen_raw(:,3);
temp = squeeze(signal_all(:,2) .* signal_all(:,3));
two_previous3_same = squeeze(signal_all((temp==1 & tr==-1)',3) .* signal_all((temp==1 & tr==-1)',4));
two_previous3_same(two_previous3_same==-1) = 0;
two_previous(2,1,1) = mean(two_previous3_same);
two_previous(2,1,2) = sqrt(var(two_previous3_same)/length(two_previous3_same));

two_previous3_diff = squeeze(signal_all((temp==-1 & tr==-1)',3) .* signal_all((temp==-1 & tr==-1)',4));
two_previous3_diff(two_previous3_diff==-1) = 0;
two_previous(2,2,1) = mean(two_previous3_diff);
two_previous(2,2,2) = sqrt(var(two_previous3_diff)/length(two_previous3_diff));

two_previous3_all = [two_previous3_same(:); two_previous3_diff(:)];
two_previous(2,3,1) = mean(two_previous3_all);
two_previous(2,3,2) = sqrt(var(two_previous3_all)/length(two_previous3_all));

two_previous_all_same = [two_previous2_same(:); two_previous3_same(:)];
two_previous(3,1,1) = mean(two_previous_all_same);
two_previous(3,1,2) = sqrt(var(two_previous_all_same)/length(two_previous_all_same));
two_previous_all_diff = [two_previous2_diff(:); two_previous3_diff(:)];
two_previous(3,2,1) = mean(two_previous_all_diff);
two_previous(3,2,2) = sqrt(var(two_previous_all_diff)/length(two_previous_all_diff));
two_previous_all = [two_previous2_all(:); two_previous3_all(:)];
two_previous(3,3,1) = mean(two_previous_all);
two_previous(3,3,2) = sqrt(var(two_previous_all)/length(two_previous_all));

% three previous case
temp = [];
tr = [];
tr = signal_chosen_raw(:,4) .* signal_not_chosen_raw(:,3);
temp = abs(squeeze(signal_all(:,1) + signal_all(:,2) + signal_all(:,3)));
three_previous3_same = squeeze(signal_all((temp==3 & tr==-1)',3) .* signal_all((temp==3 & tr==-1)',4));
three_previous3_same(three_previous3_same==-1) = 0;
three_previous(1,1) = mean(three_previous3_same);
three_previous(1,2) = sqrt(var(three_previous3_same)/length(three_previous3_same));

three_previous3_diff = squeeze(signal_all((temp==1 & tr==-1)',3) .* signal_all((temp==1 & tr==-1)',4));
three_previous3_diff(three_previous3_diff==-1) = 0;
three_previous(2,1) = mean(three_previous3_diff);
three_previous(2,2) = sqrt(var(three_previous3_diff)/length(three_previous3_diff));

three_previous3_all = [three_previous3_same(:); three_previous3_diff(:)];
three_previous(3,1) = mean(three_previous3_all);
three_previous(3,2) = sqrt(var(three_previous3_all)/length(three_previous3_all));

end