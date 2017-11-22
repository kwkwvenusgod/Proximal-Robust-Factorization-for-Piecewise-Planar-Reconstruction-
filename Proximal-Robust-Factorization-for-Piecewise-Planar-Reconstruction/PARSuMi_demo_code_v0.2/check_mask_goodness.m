function out=check_mask_goodness(mask,E,r)

mask=mask.*(~logical(E));


bad1 = sum(mask,1)<r;
bad2 = sum(mask,2)<r;

out = ~(sum(bad1)| sum(bad2));

end