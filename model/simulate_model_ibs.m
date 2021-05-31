function resp=simulate_model_ibs(theta,design_matrix)



scale_normalize=120;
orientation_std_exp=0.11;

p_match=0.7;

is_sampling=[1,1];
params=[100,theta(1),theta(2),theta(3),p_match,0,0.5,0,0.5];


%%


frame_categories=design_matrix;
frame_signals=frame_categories+normrnd(0,orientation_std_exp,size(frame_categories));
frame_signals=frame_signals*scale_normalize;



[chosen_locs,not_chosen_locs,final_choice,lo]=simulate_model(frame_signals,params,scale_normalize,is_sampling,orientation_std_exp);

resp=[chosen_locs,final_choice];

end

