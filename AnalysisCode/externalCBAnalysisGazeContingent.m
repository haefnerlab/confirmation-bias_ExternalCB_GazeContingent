function [belief_opp_cases, trials_opp] = externalCBAnalysisGazeContingent(wtchosen_signals,wtnot_chosen_signals,chosen_signals,not_chosen_signals,bias,gaps,num_peri)
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
belief_opp_all = [];
match_opp_all = [];
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
            belief_opp_all = [belief_opp_all info_per_saccade(i,j)];
            %             trials_opp = trials_opp + 1;
            match_opp_all = [match_opp_all match(i,j)];
        end
        acc_belief = acc_belief + sum(squeeze(wtnot_chosen_signals(i,:,j)));
    end
end

% trials_opp = size(belief_opp_all(:),1);
belief_opp_all = abs(belief_opp_all);
[belief_opp_all, indx_opp] = sort(belief_opp_all);
match_opp_all = match_opp_all(indx_opp);
lims = linspace(belief_opp_all(1),belief_opp_all(end),gaps);

for ii=1:size(lims,2)-1
    ind_opp_temp = find(belief_opp_all>=lims(ii) & belief_opp_all<=lims(ii+1));
    belief_opp_mid(ii) = (lims(ii)+lims(ii+1))/2;%sum(signal_opp(ind_opp_temp))/length(ind_opp_temp);
    prob_chose_in_favor_opp(ii) = mean(match_opp_all(ind_opp_temp));
    lo_err_opp(ii) = sqrt(var(match_opp_all(ind_opp_temp))/length(ind_opp_temp));
    hi_err_opp(ii) = sqrt(var(match_opp_all(ind_opp_temp))/length(ind_opp_temp));
    num_elems(ii) = length(ind_opp_temp);
end

belief_opp_cases = [belief_opp_mid; prob_chose_in_favor_opp; lo_err_opp; hi_err_opp];
trials_opp = sum(num_elems);
end