%% Script to Plot phase diagram


missing_map=ones(yn,1)*x;
outlier_map=y'*ones(1,xn);
p_map=floor((1-outlier_map).*(1-missing_map)*(m*n));
%calculate the oracle bound
O_map=OracleBound(sigma,m,n,r,p_map);



%save('LM_GN_L0_Phase_ReallyHigh2.mat')





%%

Z1_avg=median(Z1,3);
Z2_avg=median(Z2,3);
%Z3_avg=min(Z3,[],3);
Z3_avg=median(Z3,3);
%Z3_avg=Z_diff_avg;
Z4_avg=median(Z4,3);
hfig=figure(1);
set(hfig,'position',[0,0,600,500])
%image(x,y,Z3_avg,'CDataMapping','scaled')
%colormap(gray(255))
%caxis([0 10])

%Z3_avg(Z3_avg>=10)=10;

RMSE=FroNorm2RMSE(Z1_avg,m,n);
maxRMSE=FroNorm2RMSE(10,m,n);
overOracle=RMSE-O_map;
set(gca,'Position',[.09 .03 .87 .88])


image(x*100,y*100,(RMSE-O_map),'CDataMapping','scaled');
colormap(gray(255))
caxis([0 maxRMSE/2])
%caxis auto;
colorbar;
xlabel('Missing data in %','fontsize',14);
ylabel('Outliers in %','fontsize',14)
set(gca,'XAxisLocation','top','fontsize',12)
%(1-x)*m*n/dof
%axisH=axes;
%addTopXAxis(gca,'expression', '(1-argu)*m*n/dof', 'xLabStr', 'samples per d.o.f')
%axisH=axes;
%
%title('Phase diagram for recovered RMSE to ground truth')

%%




err_num_map=squeeze(sum(sum(E_gn~=0,1),2));
err_num_map=mean(err_num_map,3)+1;

E_diff_nuc=sum(sum(abs(E_gn-E_nuc3),1),2);
E_diff_avg=squeeze(mean(E_diff_nuc,5));
E_diff_avg=E_diff_avg./err_num_map;
image(E_diff_avg*50)
title('goodness of E_n_u_c initialization')
Z_diff_nuc=sqrt(sum(sum((W_gn-W_nuc3).^2,1),2));
Z_diff_avg=squeeze(mean(Z_diff_nuc,5));
figure;
image(Z_diff_avg*10)
title('goodness of W_n_u_c initialization')