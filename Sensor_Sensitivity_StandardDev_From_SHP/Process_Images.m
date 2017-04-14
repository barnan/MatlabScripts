%% startup:
clc;
clear all;
close all;
tic;

%% init:
inputFile = 'd:\MLscript\Sensor_Sensitivity_Standard_Dev_From_SHP\STD_Stack.tif';

calcFrom = 210;
calcTo = 1700;

belt1 = [500, 650];
belt2 = [1370, 1520];

%% read image:
image = imread(inputFile);

imgHeight = size(image,1);
imgWidth = size(image,2);

%figure(1), imshow(image, [min(image(:)) max(image(:))]);

%% calc bkg data:
rect1 = [10 10 50 imgHeight-20];
rect2 = [imgWidth-60 10 50 imgHeight-20];
%rect3 = [xmin ymin width height];
%rect4 = [xmin ymin width height];

imageRect1 = imcrop(image, rect1);
imageRect2 = imcrop(image, rect2);

%figure(2), imshow(imageRect1, [min(imageRect1(:)) max(imageRect1(:))]);
%figure(3), imshow(imageRect2, [min(imageRect2(:)) max(imageRect2(:))]);

bkgMean = mean2(imageRect1) + mean2(imageRect2);
bkgStd = std2(imageRect1) + std2(imageRect2);

bkgThreshold = bkgMean + 4* bkgStd;

%% calc pointlist:
pointList = [calcFrom:belt1(1,1), belt1(1,2):belt2(1,1), belt2(1,2):calcTo];

TopList1 = zeros(1, length(pointList));
TopList2 = zeros(1, length(pointList));

BotList1 = zeros(1, length(pointList));
BotList2 = zeros(1, length(pointList));

widthTop = zeros(1, length(pointList));
widthBot = zeros(1, length(pointList));

counter = 0;
for i = pointList
   counter = counter +1;
   for j = 60 : 1 : 200
       if image(j,i) > bkgThreshold
           TopList1(1,counter) = j;
           break;
       end
   end
end

counter = 0;
for i = pointList
   counter = counter +1;
   for j = 200: -1 : 60
       if image(j,i) > bkgThreshold
           TopList2(1,counter) = j;
           break;
       end
   end
end

counter = 0;
for i = pointList
   counter = counter +1;
   for j = imgHeight-60 : -1 : imgHeight-200
       if image(j,i) > bkgThreshold
           BotList1(1,counter) = j;
           break;
       end
   end
end

counter = 0;
for i = pointList
   counter = counter +1;
   for j = imgHeight-200 : 1 : imgHeight-60
       if image(j,i) > bkgThreshold
           BotList2(1,counter) = j;
           break;
       end
   end
end

%examp image:
BW = im2bw(image ./ max(image(:)), bkgThreshold ./ max(image(:)));
%figure(1); imshow(BW);

%%
widthTop = TopList2 - TopList1;
widthBot = BotList1 - BotList2;

%% find minimum:
minTop = min(widthTop(:));
minBot = min(widthBot(:));

%% position of the minimum:
posTopList = find(widthTop == minTop);
posBotList = find(widthBot == minBot);

pointTopList = pointList(posTopList);
pointBotList = pointList(posBotList);

meanPointTop = mean(pointTopList);
meanPointBot = mean(pointBotList);

%%
inputFolder = 'd:\TEMP_SUPPORT\Week03\Sensor_Sensitivity_Standard_Dev_From_SHP\';
extension = 'png';

fileList = dir([inputFolder,'*.',extension]);
kep = zeros(2048,2048);
kernel = [ 1 ; 2 ; 0 ; -2 ; -1 ];

SubPixResultList = zeros(length(fileList), 2);

for i= 1: 1 :length(fileList)
    disp([num2str(i),' / ', num2str(length(fileList))]);
    
    kep = imread([inputFolder,fileList(i).name]);
    
    kepTop = double(imcrop(kep, [meanPointTop 1 0 200 ]));
    kepBot = double(imcrop(kep, [meanPointBot imgHeight-201 0 200 ]));
    
    kepTopFiltered = imfilter(kepTop, kernel, 'replicate');
    kepBotFiltered = imfilter(kepBot, -kernel, 'replicate');
    
    maxPosTop = find(kepTopFiltered == max(kepTopFiltered(:)),1);
    maxPosBot = find(kepBotFiltered == max(kepBotFiltered(:)),1);
    
    subpixMaxPosTop = ExtrPos3(kepTopFiltered(maxPosTop-1), kepTopFiltered(maxPosTop), kepTopFiltered(maxPosTop+1), maxPosTop);
    subpixMaxPosBot = ExtrPos3(kepTopFiltered(maxPosBot+1), kepTopFiltered(maxPosBot), kepTopFiltered(maxPosBot-1), maxPosBot);
    
    SubPixResultList(i,1) = subpixMaxPosTop;
    SubPixResultList(i,2) = subpixMaxPosBot;
        
    %imshow(kepTopFiltered, [min(kepTopFiltered(:)) max(kepTopFiltered(:))]);
    %imshow(kepBotFiltered, [min(kepBotFiltered(:)) max(kepBotFiltered(:))]);
end

DeviationMean = mean(SubPixResultList);
DeviationStd = std(SubPixResultList);
DeviationRange = max(SubPixResultList) - min(SubPixResultList);

%% finish:
disp(['Elapsed time: ', num2str(toc)]);

