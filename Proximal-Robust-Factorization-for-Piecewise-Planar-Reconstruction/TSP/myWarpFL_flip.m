function [wim, count] = myWarpFL_flip(oim, flow)

% flip the flow from souce to sink

oim = double(oim);

[X, Y, C] = size(oim);

count = zeros(X,Y);
wim = zeros(X,Y,C);

for i=1:X*Y
    [x, y] = ind2sub([X Y], i);
    new_x = round(x + flow(x,y,2));
    new_y = round(y + flow(x,y,1));
    
    if (new_x > 0 && new_x <= X && new_y>0 && new_y<=Y)
        count(new_x, new_y) = count(new_x, new_y) + 1;
        wim(new_x, new_y, :) = wim(new_x, new_y, :) + oim(x, y, :);
    end
end

wim = reshape(bsxfun(@rdivide, reshape(wim, [X*Y, C]), count(:)), [X,Y,C]);