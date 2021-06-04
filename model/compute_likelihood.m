category_trial=1;
num_peri = 2;
frames = 4;
num_images = num_peri * frames +1;
num_trials=8000;

log_nsamps=linspace(0,log10(200),16);
log_sig_peri=linspace(-3,1,16);
like_x=zeros(numel(log_nsamps),numel(log_sig_peri),2^(num_images),16);
like_n=zeros(numel(log_nsamps),numel(log_sig_peri),2^(num_images),16);
for a=1:numel(log_nsamps)
    for b=1:numel(log_sig_peri)
        for c=1:2^(num_images)
            [a,b,c]
            frame_categories=(dec2bin(c-1,num_images)-'0');
            [chosen_locs,not_chosen_locs,final_choice,lo]=simulate_model_v4(repmat(frame_categories,1000,1),[100,10^log_nsamps(a),0.1,10^log_sig_peri(b),0.7,0,0.5,0,0.5],1,[1,1],0.11);
            resps=bin2dec(num2str([chosen_locs-1,0.5*(final_choice+1)]));
            [q1,q2,q3]=histcounts(resps,[0:16]-0.5);
            like_x(a,b,c,:)=q1;
            like_n(a,b,c,:)=sum(q1);
        end
    end
end