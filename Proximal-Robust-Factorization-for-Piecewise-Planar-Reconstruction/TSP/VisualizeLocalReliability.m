function  reliability=VisualizeLocalReliability(uv)
U=uv(:,:,1);
V=uv(:,:,2);
[H,W]=size(U);
reliabilityTmp=zeros(H,W);

for row=1:H
    for col=1:W
        if row<H&&row>1
            temp1=sqrt((U(row+1,col)-U(row,col))^2+(V(row+1,col)-V(row,col))^2);
            temp2=sqrt((U(row-1,col)-U(row,col))^2+(V(row-1,col)-V(row,col))^2);
            reliabilityTmp(row,col)=reliabilityTmp(row,col)+temp1+temp2;
        end
        if row==1
            temp1=sqrt((U(row+1,col)-U(row,col))^2+(V(row+1,col)-V(row,col))^2);
            reliabilityTmp(row,col)=reliabilityTmp(row,col)+temp1;
        end
        if row==H
            temp2=sqrt((U(row-1,col)-U(row,col))^2+(V(row-1,col)-V(row,col))^2);
            reliabilityTmp(row,col)=reliabilityTmp(row,col)+temp2;
        end
        if col<W&&col>1
            temp1=sqrt((U(row,col+1)-U(row,col))^2+(V(row,col+1)-V(row,col))^2);
            temp2=sqrt((U(row,col-1)-U(row,col))^2+(V(row,col-1)-V(row,col))^2);
            reliabilityTmp(row,col)=reliabilityTmp(row,col)+temp1+temp2;
        end
        if col==1
            temp1=sqrt((U(row,col+1)-U(row,col))^2+(V(row,col+1)-V(row,col))^2);
            reliabilityTmp(row,col)=reliabilityTmp(row,col)+temp1;
        end
        if col==W
             temp2=sqrt((U(row,col-1)-U(row,col))^2+(V(row,col-1)-V(row,col))^2);
             reliabilityTmp(row,col)=reliabilityTmp(row,col)+temp2;
        end
    end
    
end
 reliability=exp(-reliabilityTmp);
 
 end