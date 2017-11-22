K =180;
data='mirror';
root = ['E:\KouWen\reconstruction scene cons\data_',data,'_resize\'];
files = dir([root '*.ppm']);
dispOn = true;

[sp_labels] = TSP(K, root, files, dispOn);
savename=['E:\KouWen\reconstruction scene cons\data_',data,'_resize\','labels_',data,'_resize_K',num2str(K)];
save(savename,'sp_labels');
