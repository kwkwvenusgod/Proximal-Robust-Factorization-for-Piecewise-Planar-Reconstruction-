function check_IMG(IMG)



for k=1:IMG.K
    if (sum(IMG.label(:) == k-1) ~= IMG.SP(k).N)
        disp(['Number of pixels wrong for SP(' num2str(k) ')']);
        disp([' SP(k).N=' num2str(IMG.SP(k).N)]);
        disp([' |IMG.label==k|=' num2str(sum(IMG.label(:)==k-1))]);
        error(1);
    end
    
    [cc, Ncc] = bwlabel(IMG.label==k-1,4);
    if (Ncc>1)
        disp(['Number of connected components wrong for SP(' num2str(k) ')']);
        imagesc(cc);
        error(1);
    end
end