% Represents the loaded classifier.csv and represents the measured points from the measuredData.csv
%
% folderOfClassif   - folder of the classifier file
% nameOfClassif     - name of the represented classifier
% classesOfClassif  - the ID of the represented classes from the classifier 1-NoDefect, 2-Crack, 4-Inclusion
% colorsOfClassif   - the colors of the selected class (the order is the base of pairing
% dataOfClassif     - The represented features
%
% folderOfMeas   - folder of the classifier file
% nameOfMeas     - name of the represented classifier
% classesOfMeas  - the ID of the represented classes from the classifier 1-NoDefect, 2-Crack, 4-Inclusion
% colorsOfMeas   - the colors of the selected class (the order is the base of pairing
% dataOfMeas     - The represented features
%
% Barna N 2016
% 
%% initialization:
clear all;
close all;
tic;
clc;

%% Define input variables:

ClassifData.folderOfClassif = '.';
ClassifData.nameOfClassif = 'SemilabRuleClassifier08.csv';          % name of the classifier
ClassifData.classesOfClassif =       [ 1   2   4];                  % ID of the represented defects from the CLASSIFIER
ClassifData.colorsOfClassif =        ['r','g','b','c','y','m'];     % colors of the defects
ClassifData.dataOfClassif =          [ 0   1   2];                  % 0-DarkestGrayVal  1-MaxFilterResponse  2-Compactness

MeasData.folderOfMeas = '.';                                % path of the meas data
MeasData.nameOfMeas = 'Result_2016_2_4_11_16.csv';          % name of the Meas Data
MeasData.classesOfMeas =       [ 2 4 ];                     % ID of the represented defect from measurement
MeasData.colorsOfMeas =        ['k','m','g' ];              % colors of the defekts

Control.ShowMeas = true;

%% Read in the classifier data:

Control.runClassif = true;
ClassifData.pathOfClassif = [ClassifData.folderOfClassif,'\',ClassifData.nameOfClassif];
try
    [ClassifData.Data, ClassifData.Header] = ReadMATCSV (ClassifData.pathOfClassif);
catch exc
    Control.runClassif = false;
    msg = ['Problem during classifier file reading: ',ClassifData.pathOfClassif,' '];
    warning(msg);
end

%% read in the measeuement data:
Control.runMeas = true;
MeasData.pathOfMeas = [MeasData.folderOfMeas,'\',MeasData.nameOfMeas];

try
    [MeasData.Data, MeasData.Header] = ReadMATCSV (MeasData.pathOfMeas);
catch exc
    Control.runMeas = false;
    msg = ['Problem during measurement file reading: ',MeasData.pathOfMeas,' '];
    warning(msg);
end

%% Run or Error Message??

if ~Control.runMeas && ~Control.runClassif
    msgbox({'No such input files: ';ClassifData.pathOfClassif;MeasData.pathOfMeas}, 'Error','error');
end

%% Convert the cell arrays to double matrix

if Control.runClassif
    MatrC = zeros(size(ClassifData.Data,1), size(ClassifData.Data,2));
    if iscell(ClassifData.Data)
        for i = 1: size(ClassifData.Data,1)
            for j = 1: size(ClassifData.Data,2)
                if ~isnan(str2double(ClassifData.Data{i,j}))
                    MatrC(i,j) = str2double(ClassifData.Data{i,j});
                else
                    MatrC(i,j) = 0;
                end
            end
        end
    else
        MatrC = ClassifData.Data;
    end
end

if Control.runMeas
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
end

%% Derive variables, settings:

fulllegendOfClassif = {  '1' 'NoDefect -Classifier' ; '2' 'Crack -Classifier'; '4' 'Inclusion -Classifier'; '5' 'Edge -Classifier';...
    '6' 'Bright -Classifier'; '7' 'Double -Classifier'; '8' 'Pinhole -Classifier'; '9' 'Broken -Classifier';...
    '11' 'Sawmark -class'};

fh = findall(0,'type','figure');
str = '';

if size(MatrC,2) == 43           % 43 -> if the data was read from a classifier
    dataOffsetC = 3;
    humanCtrlColNumC = 2;
elseif size(MatrC,2) == 54       % 54 -> if the data was read from measurement csv
    dataOffsetC = 14;
    humanCtrlColNumC = 4;
else
    return;
end

datacolC = ClassifData.dataOfClassif + dataOffsetC;        % the data is expected in this columnumber

if isempty(fh)
    fig_id = 1;
else
    fig_id = max(fh.Number) + 1;
end


%% Visualize the content of the classifer:

figure(fig_id);

for i = 1: length(ClassifData.classesOfClassif)
    disp([num2str(i),' / ', num2str(length(ClassifData.classesOfClassif))]);
    
    % find indexes for the given defect:
    eval(['rowIndexC',num2str(i),' = find (MatrC(:,2) == ClassifData.classesOfClassif(',num2str(i),'));']);
    
    % copy data into represent
    eval(['DataC',num2str(i),' = [MatrC(rowIndexC',num2str(i),',',num2str(datacolC(1)),'),MatrC(rowIndexC',num2str(i),',',num2str(datacolC(2)),'), MatrC(rowIndexC',num2str(i),',',num2str(datacolC(3)),')];']);
    
    % represent:
    eval(['scatter3(DataC',num2str(i),'(:,1), DataC',num2str(i),'(:,2), DataC',num2str(i),'(:,3),ClassifData.colorsOfClassif(',num2str(i),'),''filled'' );']);
    hold all;
    
    % create legend:
    ind = find(ismember(fulllegendOfClassif,num2str(ClassifData.classesOfClassif(i))));
    str = strcat(str,'''',char(cellstr(fulllegendOfClassif{ind,2})),'''',',');
end
eval(['legend(',str(1:end-1),');']);

% if (exist('Header_selected','var') == 1)
%     head = Header_selected;
% else
%     head = Header;
% end

% label the axes:
xlabel([num2str(ClassifData.dataOfClassif(1)),' - ',char(ClassifData.Header{1,datacolC(1)})]);
ylabel([num2str(ClassifData.dataOfClassif(2)),' - ',char(ClassifData.Header{1,datacolC(2)})]);
zlabel([num2str(ClassifData.dataOfClassif(3)),' - ',char(ClassifData.Header{1,datacolC(3)})]);


%% Add the selected measurement data points to the graph:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~Control.ShowMeas
    return;
end

fulllegendMeas = {  '1' 'NoDefect -meas' ; '2' 'Crack -meas'; '4' 'Inclusion -meas'; '5' 'Edge -meas';...
    '6' 'Bright -meas'; '7' 'Double -meas'; '8' 'Pinhole -meas'; '9' 'Broken -meas';...
    '11' 'Sawmark -meas'};

% if exist ([pathOfMeas,'\',nameOfMeas,'.mat'],'file')
%     eval (['load ' pathOfMeas,'\',nameOfMeas,'.mat;']);
% elseif exist ([nameOfMeas,'.mat'],'file')
%     eval (['load ' nameOfMeas,'.mat;']);
% end

if size(MatrM,2) == 43           % 43 -> if the data was read from a classifier
    dataOffsetM = 3;
    humanCtrlColNumM = 2;
elseif size(MatrM,2) == 54       % 54 -> if the data was read from measurement csv
    dataOffsetM = 14;
    humanCtrlColNumM = 4;
else
    return;
end

datacolM = ClassifData.dataOfClassif + dataOffsetM;        % the data is expected in this columnumber

% Matr = [];
% eval(['Matr = ', nameOfMeas, ';']);

for i = 1: length(MeasData.classesOfMeas)
    disp([num2str(i),' / ', num2str(length(MeasData.classesOfMeas))]);
    
    % find indexes for the given defect:
    eval(['rowIndexM',num2str(i),' = find (MatrM(:,2) == MeasData.classesOfMeas(',num2str(i),'));']);
    
    % copy data into represent
    eval(['DataM',num2str(i),' = [MatrM(rowIndexM',num2str(i),',',num2str(datacolM(1)),'),MatrM(rowIndexM',num2str(i),',',num2str(datacolM(2)),'), MatrM(rowIndexM',num2str(i),',',num2str(datacolM(3)),')];']);
    
    % represent:
    eval(['scatter3(DataM',num2str(i),'(:,1), DataM',num2str(i),'(:,2), DataM',num2str(i),'(:,3),MeasData.colorsOfMeas(',num2str(i),'),''filled'' );']);
    hold all;
    
    % create legend:
    ind=find(ismember(fulllegendMeas,num2str(MeasData.classesOfMeas(i))));
    str = strcat(str,'''',char(cellstr(fulllegendMeas{ind,2})),'''',',');
end
eval(['legend(',str(1:end-1),');']);


%% Display of runtime
disp (['Elapsed time: ' ,num2str(toc)]);



