% Read and represent csv file features in plot:
%
% nameOfMeas:       - name of the input csv file
% showErrorBar:     - show the red errorbars or just the defect points
% 
% Barna N 2016
%
%% initialization:

clear all;
close all;
tic;
clc;

%% Define input variables:

showErrorBar = true;

MeasData.folderOfMeas = '.';                            % path of the meas data
MeasData.nameOfMeas = 'SemilabRuleClassifier08.csv';      %name of the Meas Data

%% read in the measeuement data:
MeasData.pathOfMeas = [MeasData.folderOfMeas,'\',MeasData.nameOfMeas];

try
    [MeasData.Data, MeasData.Header] = ReadMATCSV (MeasData.pathOfMeas);
catch exc
    Control.runMeas = false;
    msg = ['Problem during measurement file reading: ',MeasData.pathOfMeas,' '];
    warning(msg);
end

%%

MatrM = zeros(size(MeasData.Data,1), size(MeasData.Data,2));
if iscell(MeasData.Data)
    for i = 1: size(MeasData.Data,1)
        for j = 1: size(MeasData.Data,2)
            if isempty(MeasData.Data{i,j})
                MeasData.Data{i,j} = '0';
            end
            %                 if isnumeric(MeasData.Data{i,j});
            MatrM(i,j) = str2double(MeasData.Data{i,j});
            %                 end
        end
    end
else
    MatrM = MeasData.Data;
end


if size(MatrM,2) == 43           % 43 -> if the data was read from a classifier
    dataOffsetM = 3;
    MCIResuColumn = 2;
elseif size(MatrM,2) == 54       % 54 -> if the data was read from measurement csv
    dataOffsetM = 14;
    MCIResuColumn = 4;
else
    return;
end

%% load data:
% load data_Fortix_good_test;
% load data_MCI_Crack_0826;
% load SemilabRuleClassifier10\SemilabRuleClassifier10.mat
% load SelectedDataFromDAQ.mat
% data = data_MCI_Crack_0826;
% data = SelectedDataFromDAQ;
data_all = MatrM;
% data_ref = SemilabRuleClassifier10;

%% representation:
meas = 0; % 0 1
ReprRule = 0;

% if meas == 0
    
    %%%%% CLASSIFIER
    % Rule - 1 pre..2.3.4 post , Name - 2 crack, 4 inclu,  Featurwes - 0..40
    
    % x = rule( data_rept0903(:,1)==1, 2:end); % preclass
    % x = rule( rule(:,1)==2, 2:end);  % postclass
    % x = rule( rule(:,1)==3, 2:end);  % postclass
    
    % x = multi( multi(:,1)==1, 2:end); % preclass;
    % x = mono( mono(:,1)==1, 2:end); % preclass;
    % x = lastbasler( lastbasler(:,1)==1, 2:end); % preclass; lastbasler
    % x = bigbasler( bigbasler(:,1)==1, 2:end); % preclass; bigbasler dirt bright sawmark???
    % x = bigbasler( : , 2:end); % preclass; bigbasler dirt bright sawmark???
    
%     unk = data( data(:,2)==0,:);
%     nod = data( data(:,2)==1,:);
%     crk = data( data(:,2)==2,:);
%     inc = data( data(:,2)==4,:);
%     edg = data( data(:,2)==5,:);
%     brg = data( data(:,2)==6,:);
%     dkg = data( data(:,2)==7,:);
%     hol = data( data(:,2)==8,:);
%     drt = data( data(:,2)==9,:);
%     clu = data( data(:,2)==10,:);
%     saw = data( data(:,2)==11,:);
    
%     if ReprRule
%         nod_ref = data_ref( data_ref(:,2) == 1,:);
%         crk_ref = data_ref( data_ref(:,2) == 2,:);
%         inc_ref = data_ref( data_ref(:,2) == 4,:);
%         brg_ref = data_ref( data_ref(:,2) == 6,:);
%         dkg_ref = data_ref( data_ref(:,2) == 7,:);
%         drt_ref = data_ref( data_ref(:,2) == 9,:);
%     end
    
%     pix=30;
%     x_axis= '0';
%     y_axis= '1';
%     z_axis= '2';
%     
%     X = str2double(x_axis)+3;
%     Y = str2double(y_axis)+3;
%     Z = str2double(z_axis)+3; % in class
%     
%     figure
%     
%     %   scatter3(unk(:,X),unk(:,Y),unk(:,Z),pix,'black','filled','Marker','o');
%     %   hold on ; %figure(gcf);
%     scatter3(nod(:,X),nod(:,Y),nod(:,Z),pix,'green','filled','Marker','o');
%     hold on ; %figure(gcf);
%     scatter3(crk(:,X),crk(:,Y),crk(:,Z),pix,'red','filled','Marker','o');
%     hold on ; %figure(gcf);
%     scatter3(inc(:,X),inc(:,Y),inc(:,Z),pix,'blue','filled','Marker','o');
%     hold on ; %figure(gcf);
%     scatter3(clu(:,X),clu(:,Y),clu(:,Z),pix,'cyan','filled','Marker','o');
%     hold on ; %figure(gcf);
%     scatter3(edg(:,X),edg(:,Y),edg(:,Z),pix,'p','filled','Marker','o');
%     hold on ; %figure(gcf);
%     scatter3(saw(:,X),saw(:,Y),saw(:,Z),pix,'black','filled','Marker','o');
%     hold on ; %figure(gcf);
%     scatter3(hol(:,X),hol(:,Y),hol(:,Z),pix,'magenta','filled','Marker','o');
%     hold on ; %figure(gcf);
%     
%     legend('nod','crk','inc','clu','edg','saw','hol');
%     xlabel(x_axis);
%     ylabel(y_axis);
%     zlabel(z_axis);
%     
%     if ReprRule
%         scatter3(nod_ref(:,X),nod_ref(:,Y),nod_ref(:,Z),pix,'green','Marker','+');
%         hold on ; %figure(gcf);
%         scatter3(crk_ref(:,X),crk_ref(:,Y),crk_ref(:,Z),pix,'red','Marker','+');
%         hold on ; %figure(gcf);
%         scatter3(inc_ref(:,X),inc_ref(:,Y),inc_ref(:,Z),pix,'blue','Marker','+');
%         hold on ; %figure(gcf);
%         scatter3(dkg_ref(:,X),dkg_ref(:,Y),dkg_ref(:,Z),pix,'cyan','Marker','+');
%         hold on ; %figure(gcf);
%         scatter3(brg_ref(:,X),brg_ref(:,Y),brg_ref(:,Z),pix,'magenta','Marker','+');
%         hold on ; %figure(gcf);
%         scatter3(drt_ref(:,X),drt_ref(:,Y),drt_ref(:,Z),pix,'p','Marker','+');
%         
%         legend('nod','crk','inc','clu','edg','nod ref','crk ref','inc ref','dkg ref','brg ref','drt ref');
%     end
%     
%     hold off;
    
    
% else
%     %nod = nodef_1;
%     crk = crack_2; % meas
%     inc = incl_4;
%     
%     X=1;Y=2;Z=3; % in measure
%     
%     figure
%     % scatter3(nod(:,X),nod(:,Y),nod(:,Z),pix,'r','filled');
%     % hold on ; %figure(gcf);
%     err=24; %29;
%     scatter3(crk(1:end-err-1,X),crk(1:end-err-1,Y),crk(1:end-err-1,Z),pix,'g','filled');
%     hold on ; %figure(gcf);
%     scatter3(crk(end-err:end,X),crk(end-err:end,Y),crk(end-err:end,Z),pix,'r','filled');
%     hold on ; %figure(gcf); 24 incl error
%     title('crack');
%     hold off;
%     
%     X=1;Y=2;Z=3; % in measure
%     
%     
%     figure
%     % scatter3(nod(:,X),nod(:,Y),nod(:,Z),pix,'r','filled');
%     % hold on ; %figure(gcf);
%     
%     err=29;
%     scatter3(inc(1:end-err-1,X),inc(1:end-err-1,Y),inc(1:end-err-1,Z),pix,'g','filled');
%     hold on ; %figure(gcf);
%     scatter3(inc(end-err:end,X),inc(end-err:end,Y),inc(end-err:end,Z),pix,'r','filled');
%     hold on ; %figure(gcf); 24 incl error
%     title('incl');
%     
%     hold off;
%     
%     figure
%     % scatter3(nod(:,X),nod(:,Y),nod(:,Z),pix,'r','filled');
%     % hold on ; %figure(gcf);
%     err=24; %29;
%     scatter3(crk(1:end-err-1,X),crk(1:end-err-1,Y),crk(1:end-err-1,Z),pix,'g','filled');
%     hold on ; %figure(gcf);
%     scatter3(crk(end-err:end,X),crk(end-err:end,Y),crk(end-err:end,Z),pix,'r','filled');
%     hold on ; %figure(gcf); 24 incl error
%     
%     err=29;
%     scatter3(inc(1:end-err-1,X),inc(1:end-err-1,Y),inc(1:end-err-1,Z),pix,'b','filled');
%     hold on ; %figure(gcf);
%     scatter3(inc(end-err:end,X),inc(end-err:end,Y),inc(end-err:end,Z),pix,'m','filled');
%     hold on ; %figure(gcf); 24 incl error
%     title('crack & incl');
%     
%     hold off;
%     
% end

%% determine error bars:

for i = 1: 11
    eval(['k', num2str(i) ,' = find(data_all(:,2) == i);']);
end

% B = data_all(find(data_all(:,2) == 2),:);

for i = 3: 43
    for j = 1: 11
        eval(['AVG(i,j) = mean(data_all(k',num2str(j),',i));']);
        
        eval(['STD(i,j) = std(data_all(k',num2str(j),',i));']);
    end
end

%%
% represent:

figure;
subplot(2,2,1);  
plot (data_all(:,MCIResuColumn), data_all(:,dataOffsetM),'.'); grid on; title('0 - DarkestGrayVal');
hold all; subplot(2,2,1);  
if showErrorBar errorbar(AVG(dataOffsetM,:), STD(dataOffsetM,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

subplot(2,2,2);  
plot (data_all(:,MCIResuColumn), data_all(:,dataOffsetM +1),'.'); grid on; title('1 - MaxFilterResponse');
hold all; subplot(2,2,2);  
if showErrorBar errorbar(AVG(dataOffsetM +1,:), STD(dataOffsetM +1,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

subplot(2,2,3);  
plot (data_all(:,MCIResuColumn), data_all(:,dataOffsetM +2),'.'); grid on; title('2 - Compactness');
hold all; subplot(2,2,3);  
if showErrorBar errorbar(AVG(dataOffsetM +2,:), STD(dataOffsetM +2,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

subplot(2,2,4);  plot (data_all(:,MCIResuColumn), data_all(:,dataOffsetM +7),'.'); grid on; title('7 - LineThickNess');
hold all; subplot(2,2,4);  
if showErrorBar errorbar(AVG(dataOffsetM +7,:), STD(dataOffsetM +7,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

legend(sprintf('1 - NoDef\n2 - Crack\n4 - Inclusion\n6 - Bright\n8 - Pinhole\n10 - Cluster\n11 - Sawmark'));

figure;
subplot(2,2,1);  
plot (data_all(:,MCIResuColumn), data_all(:,11),'.'); grid on; title('8 - DarkestGrayValQuantile');
hold all; subplot(2,2,1);  
if showErrorBar errorbar(AVG(11,:), STD(11,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

subplot(2,2,2);  
plot (data_all(:,MCIResuColumn), data_all(:,12),'.'); grid on; title('9 - MaxFilterResponseQuantile');
hold all; subplot(2,2,2);  
if showErrorBar errorbar(AVG(12,:), STD(12,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

subplot(2,2,3);  
plot (data_all(:,MCIResuColumn), data_all(:,13),'.'); grid on; title('10 - MeanGrayfilteredImage');
hold all; subplot(2,2,3);  
if showErrorBar errorbar(AVG(13,:), STD(13,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

subplot(2,2,4);  plot (data_all(:,MCIResuColumn), data_all(:,14),'.'); grid on; title('11 - MedianGrayfilteredImage');
hold all; subplot(2,2,4);  
if showErrorBar errorbar(AVG(14,:), STD(14,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

legend('1 - NoDef', '2 - Crack', '4 - Inclusion', '6 - Bright', '8 - Pinhole', '10 - Cluster', '11 - Sawmark');

figure;
subplot(2,2,1);  
plot (data_all(:,MCIResuColumn), data_all(:,15),'.'); grid on; title('12 - DarkestGrayValBackground');
hold all; subplot(2,2,1);  
if showErrorBar errorbar(AVG(15,:), STD(15,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

subplot(2,2,2);  
plot (data_all(:,MCIResuColumn), data_all(:,17),'.'); grid on; title('14 - ObjectSizeValid');
hold all; subplot(2,2,2);  
if showErrorBar errorbar(AVG(17,:), STD(17,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

subplot(2,2,3);  
plot (data_all(:,MCIResuColumn), data_all(:,19),'.'); grid on; title('16 - DiagonalLenght');
hold all; subplot(2,2,3);  
if showErrorBar errorbar(AVG(19,:), STD(19,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

subplot(2,2,4);  plot (data_all(:,MCIResuColumn), data_all(:,21),'.'); grid on; title('18 - AreaNOFEATURES');
hold all; subplot(2,2,4);  
if showErrorBar errorbar(AVG(21,:), STD(21,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

legend(sprintf('1 - NoDef\n2 - Crack\n4 - Inclusion\n6 - Bright\n8 - Pinhole\n10 - Cluster\n11 - Sawmark'));

figure;
subplot(2,2,1);  
plot (data_all(:,MCIResuColumn), data_all(:,22),'.'); grid on; title('19 - WidthNOFEATURES');
hold all; subplot(2,2,1);  
if showErrorBar errorbar(AVG(22,:), STD(22,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

subplot(2,2,2);  
plot (data_all(:,MCIResuColumn), data_all(:,23),'.'); grid on; title('20 - HeightNOFEATURES');
hold all; subplot(2,2,2);  
if showErrorBar errorbar(AVG(23,:), STD(23,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

subplot(2,2,3);  
plot (data_all(:,MCIResuColumn), data_all(:,24),'.'); grid on; title('21 - BlobSizeNOFEATURES');
hold all; subplot(2,2,3);  
if showErrorBar errorbar(AVG(24,:), STD(24,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

subplot(2,2,4);  
plot (data_all(:,MCIResuColumn), data_all(:,25),'.'); grid on; title('22 - DiagonalWidth');
hold all; subplot(2,2,4);  
if showErrorBar errorbar(AVG(25,:), STD(25,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

legend(sprintf('1 - NoDef\n2 - Crack\n4 - Inclusion\n6 - Bright\n8 - Pinhole\n10 - Cluster\n11 - Sawmark'));

figure;
subplot(2,2,1);  
plot (data_all(:,MCIResuColumn), data_all(:,26),'.'); grid on; title('23 - PerimeterNOFEATURES');
hold all; subplot(2,2,1);  
if showErrorBar errorbar(AVG(26,:), STD(26,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

subplot(2,2,2);  
plot (data_all(:,MCIResuColumn), data_all(:,27),'.'); grid on; title('24 - BrightAreaNumber of Bright pixels');
hold all; subplot(2,2,2);  
if showErrorBar errorbar(AVG(27,:), STD(27,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

subplot(2,2,3);  
plot (data_all(:,MCIResuColumn), data_all(:,28),'.'); grid on; title('25 - DistanceToNearestBorder');
hold all; subplot(2,2,3);  
if showErrorBar errorbar(AVG(28,:), STD(28,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

subplot(2,2,4);  
plot (data_all(:,MCIResuColumn), data_all(:,32),'.'); grid on; title('29 - Excentricity');
hold all; subplot(2,2,4);  
if showErrorBar errorbar(AVG(32,:), STD(32,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

legend(sprintf('1 - NoDef\n2 - Crack\n4 - Inclusion\n6 - Bright\n8 - Pinhole\n10 - Cluster\n11 - Sawmark'));

figure;
subplot(2,2,1);  
plot (data_all(:,MCIResuColumn), data_all(:,33),'.'); grid on; title('30 - area contrastof stddev');
hold all; subplot(2,2,1);  
if showErrorBar errorbar(AVG(33,:), STD(33,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

subplot(2,2,2);  
plot (data_all(:,MCIResuColumn), data_all(:,34),'.'); grid on; title('31 - area contrastof mean');
hold all; subplot(2,2,2);  
if showErrorBar errorbar(AVG(34,:), STD(34,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

subplot(2,2,3);  
plot (data_all(:,MCIResuColumn), data_all(:,39),'.'); grid on; title('36 - VariationOfDarkestVals');
hold all; subplot(2,2,3);  
if showErrorBar errorbar(AVG(39,:), STD(39,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

subplot(2,2,4);  
plot (data_all(:,MCIResuColumn), data_all(:,40),'.'); grid on; title('37 - DarkDeviation');
hold all; subplot(2,2,4);  
if showErrorBar errorbar(AVG(40,:), STD(40,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

legend(sprintf('1 - NoDef\n2 - Crack\n4 - Inclusion\n6 - Bright\n8 - Pinhole\n10 - Cluster\n11 - Sawmark'));

figure;
subplot(2,2,1);  
plot (data_all(:,MCIResuColumn), data_all(:,35),'.'); grid on; title('32 - area contrastof Q1');
hold all; subplot(2,2,1);  
if showErrorBar errorbar(AVG(35,:), STD(35,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

subplot(2,2,2);  
plot (data_all(:,MCIResuColumn), data_all(:,36),'.'); grid on; title('33 - area contrastof Q2');
hold all; subplot(2,2,2);  
if showErrorBar errorbar(AVG(36,:), STD(36,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

subplot(2,2,3);  
plot (data_all(:,MCIResuColumn), data_all(:,37),'.'); grid on; title('34 - area contrastof Q3');
hold all; subplot(2,2,3);  
if showErrorBar errorbar(AVG(37,:), STD(37,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

subplot(2,2,4);  
plot (data_all(:,MCIResuColumn), data_all(:,38),'.'); grid on; title('35 - area contrastof Q4');
hold all; subplot(2,2,4);  
if showErrorBar errorbar(AVG(38,:), STD(38,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

legend(sprintf('1 - NoDef\n2 - Crack\n4 - Inclusion\n6 - Bright\n8 - Pinhole\n10 - Cluster\n11 - Sawmark'));

figure;
subplot(2,2,1);  
plot (data_all(:,MCIResuColumn), data_all(:,41),'.'); grid on; title('38 - FillFactor');
hold all; subplot(2,2,1);  
if showErrorBar errorbar(AVG(41,:), STD(41,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

subplot(2,2,2);  
plot (data_all(:,MCIResuColumn), data_all(:,42),'.'); grid on; title('39 - DimRatio');
hold all; subplot(2,2,2);  
if showErrorBar errorbar(AVG(42,:), STD(42,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

subplot(2,2,3);  
plot (data_all(:,MCIResuColumn), data_all(:,43),'.'); grid on; title('40 - DarkGrainFeature');
hold all; subplot(2,2,3);  
if showErrorBar errorbar(AVG(43,:), STD(43,:), 'xr', 'LineWidth', 2, 'MarkerSize', 14); end;

legend(sprintf('1 - NoDef\n2 - Crack\n4 - Inclusion\n6 - Bright\n8 - Pinhole\n10 - Cluster\n11 - Sawmark'));

% figure;
% subplot(2,2,1);  plot (rule(:,2), rule(:,3),'.'); grid on; title('0 - RULE - DarkestGrayVal');
% subplot(2,2,2);  plot (rule(:,2), rule(:,4),'.'); grid on; title('1 - RULE - MaxFilterResponse');
% subplot(2,2,3);  plot (rule(:,2), rule(:,5),'.'); grid on; title('2 - RULE - Compactness');




%%


