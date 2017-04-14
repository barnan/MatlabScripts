clear all;
close all;
tic;

%%
img0 = imread('d:\MLscript\WSI_Sharpness_QuickCheck\filt_0.png');
img05 = imread('d:\MLscript\WSI_Sharpness_QuickCheck\filt_05.png');
img1 = imread('d:\MLscript\WSI_Sharpness_QuickCheck\filt_1.png');
img2 = imread('d:\MLscript\WSI_Sharpness_QuickCheck\filt_2.png');
img3 = imread('d:\MLscript\WSI_Sharpness_QuickCheck\filt_3.png');
img4 = imread('d:\MLscript\WSI_Sharpness_QuickCheck\filt_4.png');
img5 = imread('d:\MLscript\WSI_Sharpness_QuickCheck\filt_5.png');

mask = [-1 -1 -1; -1 8 -1; -1 -1 -1];

%%
fimg0 = imfilter(img0, mask);
fimg05 = imfilter(img05, mask);
fimg1 = imfilter(img1, mask);
fimg2 = imfilter(img2, mask);
fimg3 = imfilter(img3, mask);
fimg4 = imfilter(img4, mask);
fimg5 = imfilter(img5, mask);

%%
val0 = mean2(fimg0(250:300,350:400));
val05 = mean2(fimg05(250:300,350:400));
val1 = mean2(fimg1(250:300,350:400));
val2 = mean2(fimg2(250:300,350:400));
val3 = mean2(fimg3(250:300,350:400));
val4 = mean2(fimg4(250:300,350:400));
val5 = mean2(fimg5(250:300,350:400));

uval0 = mean2(fimg0(260:300,1200:1250));
uval05 = mean2(fimg05(260:300,1200:1250));
uval1 = mean2(fimg1(260:300,1200:1250));
uval2 = mean2(fimg2(260:300,1200:1250));
uval3 = mean2(fimg3(260:300,1200:1250));
uval4 = mean2(fimg4(260:300,1200:1250));
uval5 = mean2(fimg5(260:300,1200:1250));




%%
figure(1); imshow(fimg1, [min(fimg1(:)) max(fimg1(:))]);
figure(2); imshow(fimg2, [min(fimg2(:)) max(fimg2(:))]);


%%
disp(['Elapsed time: ',num2str(toc)]);
