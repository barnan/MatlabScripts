%% 
clc;
clear all;
close all;
tic;

%% 
file1 = 'o:\Sorter\WSI\_UnProcessed\CPI_1_151214_FromTom\chipping\top\1214TEST1\1412201510155204QQ1214TEST1QQ2_16bit_raw.png';
file2 = 'o:\Sorter\WSI\_UnProcessed\CPI_1_151214_FromTom\chipping\top\1214TEST1\1412201510155204QQ1214TEST1QQ2.raw';

sizeX = 4096;
sizeY = 4096;

%% 
img1 = double(imread(file1));

fin = fopen(file2, 'r');
img = (fread(fin, sizeX*sizeY.*1.5, 'uint8'))';
img2 = (reshape(img, sizeY.*1.5, sizeX))';

%% 
Cimg2 = zeros(sizeX, sizeY);

for j = 1:size(Cimg2,1)
    counter = 1; 
    
    clc; disp([num2str(size(Cimg2,1)) ,'/',num2str(j)]);  
    
    for i = 1: 2 :(size(Cimg2,2))
        variab0 = img2(j,counter);
        variab1 = img2(j,counter +1);
        variab2 = img2(j,counter +2);
        counter = counter +3;
        
        valt11 = bitshift(variab0, 4);
        valt12 = bitand(variab1, 15);
        valt21 = bitshift(variab2, 4);
        valt22 = bitshift(bitand(variab1, 240), -4);
        
        v1 = bitor( valt11, valt12 );
        v2 = bitor( valt21, valt22 ); 
        
        Cimg2(j,i) = v1;
        Cimg2(j,i+1) = v2;
    end;
end;

%% subtract
resu = img1 - Cimg2;

disp(['min: ', num2str(min(resu(:)))]);
disp(['max: ', num2str(max(resu(:)))]);
disp(['mean: ', num2str(mean(resu(:)))]);
disp(['std: ', num2str(std(resu(:)))]);

%% 
figure(1); imshow(img1, [min(img1(:)) max(img1(:))]);
figure(2); imshow(img2, [min(img2(:)) max(img2(:))]);
figure(3); imshow(Cimg2, [min(Cimg2(:)) max(Cimg2(:))]);

%% 
disp(['Elapsed time: ',num2str(toc)]);
