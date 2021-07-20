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



for kk=1:numel(subjects)
kk
data_sub = readSubjectData(subjects{kk});
dmat1=data_sub.frame_categories;
dmat=dmat1(:,[1,2,5,3,6,4,7,8,9]);
resp=bin2dec(num2str([ones(numel(data_sub.choice),3),data_sub.choice(:)]));
%%
load('dump44.mat')
ll_sig_boot=get_like_sig_boot(dmat,resp,lls,1000);
sig_ll=(get_like_sig(dmat,resp,lls));

sig_lls(kk,:)=sig_ll;

[~,idd]=max(sig_ll);
sig_2s_subs(kk)=sigs(idd).^2;
hold on
plot(log10(sigs),sig_ll,'o-','linewidth',2);
drawnow;
for i=1:size(ll_sig_boot,1)
[~,idd]=max(ll_sig_boot(i,:),[],2);
sig_2_boot(kk,i)=sigs(idd).^2;
end

end