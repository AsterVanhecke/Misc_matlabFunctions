f = 'G:\seamus\development\htpalmDev\130123_ftsZHtPalm\ftsZ_htpalm_10ms_2\ftsZhtpalm_10ms_FOV15_Acq0_phPre_1\ftsZhtpalm_10ms_FOV15_Acq0_phPre_MMImages.ome.tif';

im = imread(f);
imout = preProcessPhIm(im,20,[5,5]);
imagesc(imout);

