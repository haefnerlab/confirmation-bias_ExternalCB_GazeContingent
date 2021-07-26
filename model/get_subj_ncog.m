subjects = ...
{
'GCV3-subject01'; 'GCV3-subject02';...
'GCV3-subject03'; 'GCV3-subject07';...
'GCV3-subject08'; 'GCV3-subject11';...
'GCV3-subject13'; 'GCV3-subject15';...
'GCV3-subject16'; 'GCV3-subject17';...
'GCV3-subject18'; 'GCV3-subject19';...
'GCV3-subject20'; 'GCV3-subject21';...
'GCV3-subject22'; 'GCV3-subject23';...
};



for kk=[1,2,3,8,9,10,12,13,15,16]
kk
data_sub = readSubjectData(subjects{kk});
dmat1=data_sub.frame_categories;
dmat=dmat1(:,[1,2,5,3,6,4,7,8,9]);
resp=bin2dec(num2str([ones(numel(data_sub.choice),3),data_sub.choice(:)]));
%%

id=find(sqrt(sig_2s_subs(kk))==sqrt(sig_subs));
lls_sub=lls(:,:,:,id);

ll_cog=get_like_ncog(dmat,resp,lls_sub);
subplot(4,4,kk);
plot(ns,ll_cog,'o-','linewidth',2)
drawnow;
ll_cogs(kk,:)=ll_cog;

end