
opt='';
num_para=108;

opt='_k';
num_para=8;
DD=['/csail/fisher5/projects/Vseg_sp/outputs' opt '/'];
names={'camera_motion1','sample4','sample1','sample3',...
'RubberWhale', 'Hydrangea', 'Dimetrodon', 'Grove2', 'Grove3', 'Urban2', 'Urban3','Venus'};

re_flo=zeros(numel(names),num_para,2);
for db_id = 1:num_para
        load([DD sprintf('%06d/',db_id) 'statsf_MIT'])
        for i=1:length(flo_mis)
            re_flo(i,db_id,1)=mean(flo_mis{i})/prod(osizes(i,1:2));
            re_flo(i,db_id,2)=osizes(i,3);
        end
        load([DD sprintf('%06d/',db_id) 'statsf_MID'])
        for i=1:length(flo_mis)
            re_flo(4+i,db_id,1)=mean(flo_mis{i})/prod(osizes(i,1:2));
            re_flo(4+i,db_id,2)=osizes(i,3);
        end
end

eval(['save vseg_sp_flo' opt ' re_flo'])
