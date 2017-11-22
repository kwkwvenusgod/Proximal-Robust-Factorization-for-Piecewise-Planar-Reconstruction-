function [IMG,del_f] = SP_gen(IMG,img,dispOn)
    it = 0;
    converged = false;
    del_f=[];
    
    while (~converged) 
        it = it + 1;

        oldK = IMG.K;
        IMG.SP_changed(:) = true;
        [IMG.K, IMG.label, IMG.SP, IMG.SP_changed, IMG.max_UID,del] = local_move(IMG,100);
        del_f=[del_f,del];
        converged = sum(IMG.SP_changed);

        [IMG.K, IMG.label, IMG.SP, IMG.SP_changed, IMG.max_UID,del] = merge_move(IMG,1);
        del_f(end) = del;
        converged = converged + sum(IMG.SP_changed);
        
        %pre_l = IMG.label;    
                
        [IMG.K, IMG.label, IMG.SP, IMG.SP_changed, IMG.max_UID,del] = split_move(IMG,1);
        del_f(end) = del;        
        %{%}
        
        converged = (converged + sum(IMG.SP_changed))==0;
        if (dispOn>=1)
            sfigure(1);
            subplot(1,1,1);
            %im = zeros(size(img));
            %im(:) = c(IMG.label+1,:);
            imagesc(IMG.label);
            %image(im,'parent',gca);
            title([num2str(it) ' - ' num2str(numel(unique(IMG.label)))]);
        end
         if (dispOn>=2)
            sfigure(2);
            subplot(1,1,1);
            im = double(img)/255;
            borders = is_border_valsIMPORT(double(reshape(IMG.label+1, [size(img,1), size(img,2)])));
            im = setPixelColors(im, find(borders), [0 1 0]);
            image(im,'parent',gca);
            drawnow;
            if(dispOn==3)waitforbuttonpress;end
         end
        %sfigure(3);imagesc(pre_l)
    end
