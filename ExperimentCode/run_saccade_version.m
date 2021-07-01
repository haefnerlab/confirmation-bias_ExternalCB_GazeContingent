clear all;close all;clc;
is_new = 1;
if is_new==0
    subjectID = '';
end
%% fixed for both staircase and main experiment
num_im = 3;
con = 20;
ratio = 0.7;
track_eye = 1;
%% staircase specific variables
phase = 2; % if noise staircase
perf_lo = 0.6;
perf_hi = 0.7;
%%


if is_new==1
    proto_staircase = newGaborData('trials_per_block',200,'contrast',con,'blocks',1,'number_of_images', num_im, 'no_staircase',0, 'stair_fn', @Staircase.noise, 'ratio', ratio, 'eyelink_use', track_eye);
    params_staircase = ExperimentGabor_staircase(proto_staircase);
    subjectID = params_staircase.subjectID;
    [floor, thresh, ~] = GaborAnalysis.getThresholdWindow(subjectID, phase, perf_lo, perf_hi);
else

     params_staircase = LoadAllSubjectData(subjectID,phase);
     subjectID = params_staircase.subjectID;
    [floor, thresh, ~] = GaborAnalysis.getThresholdWindow(subjectID, phase, perf_lo, perf_hi);
    proto_exp = newGaborData('trials_per_block',100,'contrast',con,'blocks',10,'number_of_images', num_im, 'no_staircase',1,'ratio', ratio, 'eyelink_use', track_eye, 'noise', thresh);
    GaborData = ExperimentGabor(proto_exp);
end