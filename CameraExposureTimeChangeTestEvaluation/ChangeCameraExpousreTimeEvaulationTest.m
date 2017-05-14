%% initialize:
clear all;
close all;
clc;
tic;

%% define input folder:
inputFolder = 'images';

%% read in used exposure times:
fileList = dir([inputFolder, '\\*.png']);
msList = zeros(length(fileList),1);

for i = 1: length(fileList)
    [ms, remain] = strsplit(fileList(i).name,'ms');
    
    msList(i) = str2double(ms{1});
end


%% read in image and pre-process (make a single line form the center of the image)

SummDark =  [];
SummBright = [];
tempD = [];
tempB = [];
counter = 1;

for m = 1: length(fileList)
    
    clc;
    disp([num2str(m), ' / ', num2str(length(fileList))]);
        
    Dark = zeros(1,60);
    Bright = zeros(1,60);
    
    kep = imread([inputFolder, '\\', fileList(m).name]);
    threshold = 100;
    
    staX = size(kep,2) / 4;
    staY = 800; %size(kep,1) / 4;
    widthX = size(kep,2) / 2 ;
    heightY = size(kep,1) -1600; %size(kep,1) / 4;
    
    kepSmall = imcrop(kep, [staX staY widthX heightY]);
    vectVertical = mean(kepSmall, 2);
    
    % edge search and segmentation:
    filter = [-1 ; 1];
    filteredVectVertical = imfilter(vectVertical, filter, 'replicate');
    
    DarkToBrightVector = filteredVectVertical;
    DarkToBrightVector(filteredVectVertical < threshold) = 0;
    
    BrightToDarkVector = filteredVectVertical;
    BrightToDarkVector(filteredVectVertical > -threshold) = 0;
    
    BrightToDark = find (BrightToDarkVector);
    DarkToBright = find (DarkToBrightVector);
    
    %count the step size:
    
    for i = 1: (min(length(BrightToDark), length(DarkToBright))-1)
        
        if (DarkToBright(1) < BrightToDark(1))
            Bright(i) = BrightToDark(i) - DarkToBright(i);
            Dark(i) = DarkToBright(i+1) - BrightToDark(i);
        else
            Bright(i) = BrightToDark(i+1) - DarkToBright(i);
            Dark(i) = DarkToBright(i) - BrightToDark(i);
        end
        
    end
    
    tempB = horzcat(tempB, Bright);
    tempD = horzcat(tempD, Dark);
    
    if ((m < length(fileList)-1) && (msList(m) ~= msList(m+1))) || (m == length(fileList))
        
        % axisXValues will contain the ms values
        axisX(counter,1) = msList(m);
        
        counter = counter + 1;
        
        SummDark = vertcat(SummDark, tempD);
        SummBright = vertcat(SummBright, tempB);
        
        tempD = [];
        tempB = [];    
    end
    
    
end

%% calculate mean and std:

meanD = zeros(size(SummDark,1), 1);
meanB = zeros(size(SummBright,1), 1);
stdD = zeros(size(SummDark,1), 1);
stdB = zeros(size(SummBright,1), 1);

for m = 1: size(SummDark,1)
    
    tempD = SummDark(m, find(SummDark(m,:)));
    tempB = SummBright(m, find(SummBright(m,:)));
    
    meanD(m) = mean(tempD,2);
    meanB(m) = mean(tempB,2);    
        
    stdD(m) = std(tempD,0,2);
    stdB(m) = std(tempB,0,2);        
end

%% represent: 
figure (1);
errorbar(axisX, meanB, stdB *3); 
title('BRIGHT Step-size vs. Delay');
xlabel('Delay of exposure change [ms] ');
ylabel('Step size [pixel]');

figure (2);
errorbar(axisX, meanD, stdD *3, 'red'); 
title('DARK Step-size vs. Delay');
xlabel('Delay of exposure change [ms] ');
ylabel('Step size [pixel]');


%% close
disp(['elapsed time:', num2str(toc)]);

