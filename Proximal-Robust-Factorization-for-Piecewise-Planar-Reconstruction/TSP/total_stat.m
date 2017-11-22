


opt='';
num_para=108;
opt='_k';

num_para=8;
fns = dir(['output_stat' opt '/*.mat']);

stat_seg=zeros(numel(fns)*5,num_para);
for i=1:numel(fns)
    load(['output_stat' opt '/' fns(i).name]);    
    stat_seg(1+5*(i-1):5*i,:) = result;
end
eval(['save vseg_sp_stat' opt ' stat_seg'])
