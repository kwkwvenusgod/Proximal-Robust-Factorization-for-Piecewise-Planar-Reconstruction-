function [ patchflow, patchflow_rank4,label_num,structure_mask,residual] = Affine_Estimation(denseflow,label_ref,patch_similarity)
    %SMOOTHNESS_ESTIMATE_FLOW_GRADIENT Summary of this function goes here
    %   Detailed explanation goes here
    % Estimate neighbour structure

    options.epsilon = 1e-6;
    options.P_inlier = 0.99;
    options.sigma = 8*1e-3;
    options.est_fun = @estimate_plane;
    options.man_fun = @error_plane;
    options.mode = 'MSAC';
    options.Ps = [];
    options.notify_iters = [];
    options.min_iters = 50;
    options.fix_seed = false;
    options.reestimate = true;
    options.stabilize = false;
    % image coordinate reset
    discardL=[];
    DISCARD=[];
    LL=[];
    labelmap=label_ref;
    [H,W]=size(label_ref);
    label_num=unique(labelmap);
    ind=find(label_num==inf);
    label_num(ind)=[];
   
    smoothness_structure=zeros(6*length(label_num));
    for i=1:length(label_num)
        tmp=patch_similarity(:,i);
        tmp=repmat(tmp,1,6);
        tmp=reshape(tmp',6*length(label_num),1);
        tmp=repmat(tmp,1,6);
        smoothness_structure(:,i*6-5:i*6)=tmp;
    end
    numofframes=size(denseflow,2);
    discard_count=0;
    discard_count_global=0;
    beta=75;

    for i=1:numofframes
        label=label_ref;
        dim_A=0;
        flow=denseflow{i};
        flowu=-flow(:,:,1);
        flowv=-flow(:,:,2);

        label_num_tmp=label_num;
        smoothness_construct_tmp=smoothness_structure;
        patchrank=ones(length(label_num),1);
        tic
        for j=1:length(label_num)

            ll=label_num(j);
            ind=find(label==ll);
            [y_uv,x_uv]=find(label==ll);
            local_u=flowu(ind);
            local_v=flowv(ind);
            x_uv=-(x_uv-W/2*ones(size(x_uv)));
            y_uv=H/2*ones(size(y_uv))-y_uv;
            Datau=[local_u';x_uv';y_uv'];
            Datav=[local_v';x_uv';y_uv'];
            %        X=[X1;X2];
            if size(Datau,2)>10
                [results, options] = RANSAC(Datau, options);
                RANSAC_mask_u=results.CS;
                %                        eu=mean(results.E);
                [results, options] = RANSAC(Datav, options);
                %                        ev=mean(results.E);

                RANSAC_mask_v=results.CS;
                RANSAC_mask=and(RANSAC_mask_u,RANSAC_mask_v);

                ind_ransac=find(RANSAC_mask==1);
                local_u_ransac=local_u(ind_ransac);
                local_v_ransac=local_v(ind_ransac);
                x_local=x_uv(ind_ransac);
                y_local=y_uv(ind_ransac);



                tmp=[ones(length(x_local),1),x_local,y_local];
                A=zeros(length(x_local)*2,6);
                A(1:length(x_local),1:3)=tmp;
                A(1+length(y_local):end,4:6)=tmp;
                value_uv=[local_u_ransac;local_v_ransac];%  eu>0.4||ev>0.4

                %              ||length(ind_ransac)/size(Datau,2)<0.7

                if rank(A)<6
                    discard_count=discard_count+1;
                    discard_count_global=discard_count_global+1;
                    discardL(discard_count)=label_num(j);
                    DISCARD(discard_count_global)=label_num(j);
                    continue;
                else
                    patchrank(j)=length(ind_ransac)/size(Datau,2);
                    %                   affinetmp=pinv(A)*value_uv;
                    %                   tmp=[ones(length(x_uv),1),x_uv,y_uv];
                    %                   A_all=zeros(length(x_uv)*2,6);
                    %                   A_all(1:length(x_uv),1:3)=tmp;
                    %                   A_all(1+length(x_uv):end,4:6)=tmp;
                    %                   b_all=[local_u;local_v];
                    %                   restmp=norm(A_all*affinetmp-b_all);
                    %                   patchrank(j)=restmp;
                end
                dim_A=dim_A+length(x_local)*2;

            else
                discard_count=discard_count+1;
                discard_count_global=discard_count_global+1;
                discardL(discard_count)=label_num(j);
                DISCARD(discard_count_global)=label_num(j);
                continue;
            end
            % tmp=[ones(length(x_uv),1),x_uv,y_uv];
            subA{j-discard_count}=A;
            subb{j-discard_count}=value_uv;
            LL(j-discard_count)=label_num(j);
        end

        %     [C,inx]=sort(patchrank,'descend');
        %     [C,inx]=sort(patchrank);
        %     for k=1:numtmp
        %         discard_count=discard_count+1;
        %         discardL(discard_count)=label_num(inx(k));
        %         discard_count_global=discard_count_global+1;
        %         DISCARD(discard_count_global+1)=label_num(inx(k));
        %     end
        [discard_value,ind_dicard]=intersect(label_num,discardL);
        label_num_tmp(ind_dicard)=[];
        A=zeros(dim_A,length(label_num_tmp)*6);
        b=zeros(dim_A,1);
        startingrow=1;
        for k=1:length(label_num_tmp)
            ind=find(LL==label_num_tmp(k));
            if isempty(ind)
                break;
            end
            tmp=subA{ind};
            n=size(tmp,1);
            A(startingrow:startingrow+n-1,k*6-5:k*6)=tmp;
            tmp=subb{ind};
            b(startingrow:startingrow+n-1)=tmp;
            startingrow=startingrow+n;
        end


        for k=1:length(ind_dicard)
            cc=ind_dicard(k)-k+1;
            smoothness_construct_tmp(cc*6-5:cc*6,:)=[];
            smoothness_construct_tmp(:,cc*6-5:cc*6)=[];
        end

        A=sparse(A);
        C=A'*A+beta*smoothness_construct_tmp;
        d=A'*b;
        parai=(C)\d;
        L{i}=label_num_tmp;
        PARA{i}=parai;
        toc
        discard_count=0;
        discardL=[];
        LL=[];
    end
    % discard corrupted patches
    DISCARD=unique(DISCARD);
    tmp=find(DISCARD==inf);
    DISCARD(tmp)=[];
    numofpatches=length(label_num)-length(DISCARD);


    patchflow=zeros(numofframes,6*numofpatches);
    patchflow_rank4=zeros(numofframes,6*numofpatches);
    for i=1:numofframes
        ll= L{i};
        tmp=PARA{i};
        tmp=tmp';
        for j=1:length(DISCARD)
            dd=DISCARD(j);
            ind=find(ll==dd);
            if isempty(ind)==0
                tmp(:,6*ind-5:6*ind)=[];
                ll(ind)=[];
            end

        end
        tmp1=reshape(tmp,6,length(tmp)/6);
        [u,s,v]=svds(tmp1',4);
        tmp1=u*s*v';
        tmp1=tmp1';
        tmp1=tmp1(:);
        patchflow_rank4(i,:)=tmp1';
        patchflow(i,:)=tmp';
        % patchflow(i,:)=tmp;
    end

    [tmp,tmp1]=intersect(label_num,DISCARD);
    label_num(tmp1)=[];
    %rearrange the matrix
    patchflow_rearrange=zeros(numofframes,6*numofpatches);
    row=1:numofpatches;
    r=1:numofframes;
    r=repmat(r,numofpatches,1);
    r=reshape(r',numofpatches*numofframes,1);


    tmp=patchflow(:,row*6-5*ones(1,numofpatches));
    tmp=tmp(:);
    c=row*6-5*ones(1,numofpatches);
    c=repmat(c,numofframes,1);
    c=reshape(c,numofframes*numofpatches,1);
    ind=sub2ind(size(patchflow_rearrange),r,c);
    patchflow_rearrange(ind)=tmp;

    tmp=patchflow(:,row*6-4*ones(1,numofpatches));
    tmp=tmp(:);
    c=row*6-3*ones(1,numofpatches);
    c=repmat(c,numofframes,1);
    c=reshape(c,numofframes*numofpatches,1);
    ind=sub2ind(size(patchflow_rearrange),r,c);
    patchflow_rearrange(ind)=tmp;

    tmp=patchflow(:,row*6-3*ones(1,numofpatches));
    tmp=tmp(:);
    c=row*6-1*ones(1,numofpatches);
    c=repmat(c,numofframes,1);
    c=reshape(c,numofframes*numofpatches,1);
    ind=sub2ind(size(patchflow_rearrange),r,c);
    patchflow_rearrange(ind)=tmp;

    tmp=patchflow(:,row*6-2*ones(1,numofpatches));
    tmp=tmp(:);
    c=row*6-4*ones(1,numofpatches);
    c=repmat(c,numofframes,1);
    c=reshape(c,numofframes*numofpatches,1);
    ind=sub2ind(size(patchflow_rearrange),r,c);
    patchflow_rearrange(ind)=tmp;

    tmp=patchflow(:,row*6-1*ones(1,numofpatches));
    tmp=tmp(:);
    c=row*6-2*ones(1,numofpatches);
    c=repmat(c,numofframes,1);
    c=reshape(c,numofframes*numofpatches,1);
    ind=sub2ind(size(patchflow_rearrange),r,c);
    patchflow_rearrange(ind)=tmp;

    tmp=patchflow(:,row*6-0*ones(1,numofpatches));
    tmp=tmp(:);
    c=row*6-0*ones(1,numofpatches);
    c=repmat(c,numofframes,1);
    c=reshape(c,numofframes*numofpatches,1);
    ind=sub2ind(size(patchflow_rearrange),r,c);
    patchflow_rearrange(ind)=tmp;

    patchflow=patchflow_rearrange;


    %rearrange the matrix
    patchflow_rearrange=zeros(numofframes,6*numofpatches);
    row=1:numofpatches;
    r=1:numofframes;
    r=repmat(r,numofpatches,1);
    r=reshape(r',numofpatches*numofframes,1);


    tmp=patchflow_rank4(:,row*6-5*ones(1,numofpatches));
    tmp=tmp(:);
    c=row*6-5*ones(1,numofpatches);
    c=repmat(c,numofframes,1);
    c=reshape(c,numofframes*numofpatches,1);
    ind=sub2ind(size(patchflow_rearrange),r,c);
    patchflow_rearrange(ind)=tmp;

    tmp=patchflow_rank4(:,row*6-4*ones(1,numofpatches));
    tmp=tmp(:);
    c=row*6-3*ones(1,numofpatches);
    c=repmat(c,numofframes,1);
    c=reshape(c,numofframes*numofpatches,1);
    ind=sub2ind(size(patchflow_rearrange),r,c);
    patchflow_rearrange(ind)=tmp;

    tmp=patchflow_rank4(:,row*6-3*ones(1,numofpatches));
    tmp=tmp(:);
    c=row*6-1*ones(1,numofpatches);
    c=repmat(c,numofframes,1);
    c=reshape(c,numofframes*numofpatches,1);
    ind=sub2ind(size(patchflow_rearrange),r,c);
    patchflow_rearrange(ind)=tmp;

    tmp=patchflow_rank4(:,row*6-2*ones(1,numofpatches));
    tmp=tmp(:);
    c=row*6-4*ones(1,numofpatches);
    c=repmat(c,numofframes,1);
    c=reshape(c,numofframes*numofpatches,1);
    ind=sub2ind(size(patchflow_rearrange),r,c);
    patchflow_rearrange(ind)=tmp;

    tmp=patchflow_rank4(:,row*6-1*ones(1,numofpatches));
    tmp=tmp(:);
    c=row*6-2*ones(1,numofpatches);
    c=repmat(c,numofframes,1);
    c=reshape(c,numofframes*numofpatches,1);
    ind=sub2ind(size(patchflow_rearrange),r,c);
    patchflow_rearrange(ind)=tmp;

    tmp=patchflow_rank4(:,row*6-0*ones(1,numofpatches));
    tmp=tmp(:);
    c=row*6-0*ones(1,numofpatches);
    c=repmat(c,numofframes,1);
    c=reshape(c,numofframes*numofpatches,1);
    ind=sub2ind(size(patchflow_rearrange),r,c);
    patchflow_rearrange(ind)=tmp;

    patchflow_rank4=patchflow_rearrange;

    mask=isnan(patchflow);
    [r,c]=find(double(mask)==1);
    discart_list=unique(c);
    if ~isempty(discart_list)
        discard_label2=fix(double((discart_list-ones(length(discart_list),1))/6))+ones(length(discart_list),1);
        discard_label2=unique(discard_label2);
        count=0;
        for i=1:length(discard_label2)
            k=discard_label2(i)-count;
            patchflow(:,k*6-5:k*6)=[];
            patchflow_rank4(:,k*6-5:k*6)=[];
            label_num(k)=[];
            count=count+1;
        end
    end

    structure_mask=inf*ones(size(labelmap));
    for i=1:length(label_num)
        ind=find(labelmap==label_num(i));
        structure_mask(ind)=label_num(i);
    end

    %calculate residual
    residual=zeros(size(patchflow,1),size(patchflow,2)/6);
    for i=1:size(patchflow,1)
        tmpuv=denseflow{i};
        tmpu=tmpuv(:,:,1);
        tmpv=tmpuv(:,:,2);
        tmplable=label_ref;
        for j=1:size(patchflow,2)/6
            ll=label_num(j);
            [y,x]=find(tmplable==ll);
            ind=find(tmplable==ll);
            uu=tmpu(ind);
            vv=tmpv(ind);
            b=[uu;vv];
            tmp=[ones(length(x),1),x,y];
            A=zeros(length(x)*2,6);
            A(1:length(x),1:3)=tmp;
            A(length(x)+1:end,4:end)=tmp;
            para=patchflow(i,j*6-5:j*6);
            para=[para(1),para(3),para(5),para(2),para(4),para(6)];
            output=A*para';
            r=output-b;
            residual(i,j)=norm(r)/length(output);
        end
    end

end

