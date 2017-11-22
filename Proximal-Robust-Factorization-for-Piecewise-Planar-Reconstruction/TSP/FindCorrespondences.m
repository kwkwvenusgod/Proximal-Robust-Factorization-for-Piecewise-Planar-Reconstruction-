function [ correspondenceMatrix  ] = FindCorrespondences( labels)
%FINDCORRESPONDENCES Summary of this function goes here
%   Detailed explanation goes here

[H, W, numofframes]=size(labels);

correspondenceMatrix=inf*ones(H, W ,numofframes);

stat=unique(labels(:,:,1));

record_d_label=1;
for count=1:length(stat)
    for f=2:numofframes
        l=labels(:,:,f);
        [row, col]=find(l==stat(count));
        if isempty(row)
            d_l(record_d_label)=stat(count);
            record_d_label=record_d_label+1;
        end
    end
end
d_ll=unique(d_l);
for i=1:length(d_ll)
    pos=find(stat==d_ll(i));
    stat(pos)=[];
end

for i=1:numofframes
    l=labels(:,:,i);
    for j=1:length(stat)
        [row, col]=find(l==stat(j));        
        ind=sub2ind(size(labels),row,col,repmat(i,length(row),1));
        correspondenceMatrix(ind)=stat(j); 
    end
end

end

