

h = imagesc(sp_labels(:,:,1),[0 14405051]);
freezeColors;
im = get(h, 'CData');
imwrite(uint8(im*255), '~/Desktop/Latex/jchang7_2013_cvpr_sp/figures/table_swa1.png');

h = imagesc(sp_labels(:,:,2),[0 14405051]);
freezeColors;
im = get(h, 'CData');
imwrite(uint8(im*255), '~/Desktop/Latex/jchang7_2013_cvpr_sp/figures/table_swa2.png');