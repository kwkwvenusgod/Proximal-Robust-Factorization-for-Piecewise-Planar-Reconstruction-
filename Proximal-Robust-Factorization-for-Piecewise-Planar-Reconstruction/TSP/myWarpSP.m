function [wim, count] = myWarpSP(oim, label, vx, vy, K)

X = size(oim,1);
Y = size(oim,2);
C = size(oim,3);

wim = zeros(size(oim));
count = zeros(X,Y);

for k=1:K
    mask = (label==k);
    
    mask_indices = find(mask);
    [xi,yi] = ind2sub([X,Y],mask_indices);
    x = round(xi + vx(k));
    y = round(yi + vy(k));
    
    valid = (x>0 & x<=X & y>0 & y<=Y);
    x = x(valid);
    y = y(valid);
    xi = xi(valid);
    yi = yi(valid);
    
    if (~isempty(x) && ~isempty(y))
        input_indices = sub2ind([X,Y], xi, yi);
        output_indices = sub2ind([X,Y], x, y);
        count(output_indices) = count(output_indices) + 1;

        offset = repmat([0:X*Y:X*Y*(C-1)], [size(output_indices,1), 1]);
        output_indices = repmat(output_indices, [1 C]);
        output_indices = output_indices + offset;
        input_indices = repmat(input_indices, [1 C]);
        input_indices = input_indices + offset;

        wim(output_indices(:)) = wim(output_indices(:)) + oim(input_indices(:));
    end
end

wim = reshape(bsxfun(@rdivide, reshape(wim, [X*Y, C]), count(:)), [X,Y,C]);
wim = setPixelColors(wim, isnan(wim(:,:,1)), [255 0 0]);