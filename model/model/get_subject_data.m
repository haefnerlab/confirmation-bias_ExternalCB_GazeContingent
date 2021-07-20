sessions=[1,2,3];
orientation=45;
num_frames=4;
num_images=2*num_frames+1;
sub_struct.num_frames=num_frames;
sub_struct.signals=[];
sub_struct.saccades=[];
sub_struct.correct_resp=[];
sub_struct.choice=[];
sub_struct.num_trials=0;
for i=sessions
    i
    load(['../RawData/GCV3-subject03-Session',num2str(i),'-GaborDataNoStaircaseQuit.mat']);
    ids=find(GaborData.total_number_of_images==num_images);
    for j=ids
        

        sub_struct.signals=[sub_struct.signals;GaborData.ideal_frame_signals{j}(:)'];
        tmp1=GaborData.image_array_chosen_index{j}(:)';
        sub_struct.saccades=[sub_struct.saccades;(1+(GaborData.image_array_chosen_index{j}(2:end-1,2)==1))'];
        sub_struct.correct_resp=[sub_struct.correct_resp;sign(GaborData.ideal_answer(j)-0.5)];
        sub_struct.choice=[sub_struct.choice;sign(GaborData.choice(j)-0.5)];
        sub_struct.num_trials=sub_struct.num_trials+1;
    end
end
