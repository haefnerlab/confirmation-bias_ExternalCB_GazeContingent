function [prob_opp_divided] = externalCBAnalysisGazeContingent_divided(wtchosen_signals,wtnot_chosen_signals,chosen_signals,not_chosen_signals,bias,num_peri)
% This function returns probability of choosing to saccade in a location
% favoring evidence collected so far
% Input: chosen_signals = array of chosen signal values, trials x number of
% frames
% not_chosen_signals = array of  not chosen signal values, trials x number of
% frames-1
% we hope these signals are appropriately weighted

% Output: belief_mid is binned values of accumulated evidence and
% chose_in_favor is corresponding probability of choosing evidence in favor
% of accumulation

num_trials = size(chosen_signals,1);
count_opp = zeros(1,size(chosen_signals,2)-1);
match_opp_all = zeros(1,size(chosen_signals,2)-1);
prob_opp_divided = zeros(1,size(chosen_signals,2)-1);
trials_opp = 0;
num_saccade = size(chosen_signals,2)-1;
for i=1:num_trials
    acc_belief = bias;
    for j=1:num_saccade % evidence guiding saccade is collected till max saccade - 1
        info_per_saccade(i,j) = acc_belief + wtchosen_signals(i,j) ;% evidence collected so far
        evidence_next = [chosen_signals(i,j+1) squeeze(not_chosen_signals(i,:,j))];
        acc_belief = info_per_saccade(i,j);
        match(i,j) = sign(evidence_next(1))==sign(acc_belief);
        
        if (abs(sum(sign(evidence_next)))~=num_peri)
            count_opp(j)  = count_opp(j) + 1;
            match_opp_all(j) = match_opp_all(j) + match(i,j);
        end
        acc_belief = acc_belief + sum(squeeze(wtnot_chosen_signals(i,:,j)));
    end
end
for j=1:num_saccade
    if count_opp(j)>0
        prob_opp_divided(j,1) = match_opp_all(j)/count_opp(j);
        prob_opp_divided(j,2) = sqrt((prob_opp_divided(j,1)*(1-prob_opp_divided(j,1))))/sqrt(count_opp(j));
    else
        prob_opp_divided(j,1) = 0;
        prob_opp_divided(j,2) = 0;
    end
end
end