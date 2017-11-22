function M_list=construct_M_list(M,mask)
%transform M into cell list form
[m,n]=size(mask);
M_list=cell(n,1);
idx=sum(mask);
idx1=1;
idx2=idx(1);
M_list{1}=M(idx1:idx2);
for i=2:n
    idx1=1+idx2;
    idx2=idx2+idx(i);
    M_list{i}=M(idx1:idx2);
end

end