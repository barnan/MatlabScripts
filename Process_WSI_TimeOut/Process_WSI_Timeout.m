
%% initialize:
clear all;
close all;
tic;
clc;

%% read csv file:
filename1 = 'd:\MLscript\ProcessWSI_Timeout\WINS_11__2016090400.csv';

disp(['file reading: ', filename1]);

[Data Text] = xlsread(filename1);

%% process the header:
rem = char(Text{1,1});
counter = 1;

Header = cell(1,408);

for i = 1 : size(Header,2)
    [val, rem] = strtok(rem,',');
    Header{1,i} = val;
end

disp ('Header converted');

% define  inspected columns:
column.PLI = 24;
column.SHP = 107;
column.MCI = 109;
column.TTR = 365;
column.ChipTop = 143;
column.ContTop = 221;
column.ChipBot = 265;
column.ContBot = 329;
column.Lost = 6;

%%  process the content of the csv file:
% Content = cell (size(Text,1)-1, size(Header,2));

counter.ChipTop = 0;
counter.ChipBot = 0;
counter.ContTop = 0;
counter.ContBot = 0;

Content = cell(1,408);

for j = 2: (size(Text,1))
    rem = Text{j,1};
    
    if mod(j,1000) == 0
        clc; disp(['row:', num2str(j)]);
    end
    
    indexes = strfind(Text{j},',');
    
    cont = Text{j,1};
    
    for i = 1 : (size(indexes,2))
        
        if i > 1
            index1 = indexes(i-1)+1;
        else
            index1 = 1;
        end
        
        index2 = indexes(i)-1;
        
        if index1 == (index2+1)
            Content{1,i} = '';
        else
            Content{1,i} = cont(index1 :index2);
        end               
    end
    
    if isempty(Content{1,column.ChipTop})
        counter.ChipTop = counter.ChipTop +1;
    end
    if isempty(Content{1,column.ChipBot})
        counter.ChipBot = counter.ChipBot +1;
    end
    if isempty(Content{1,column.ContTop})
        counter.ContTop = counter.ContTop +1;
    end
    if isempty(Content{1,column.ContBot})
        counter.ContBot = counter.ContBot +1;
    end       
end

%% Study WSI columns
% counter.ChipTop = 0;
% counter.ChipBot = 0;
% counter.ContTop = 0;
% counter.ContBot = 0;
% % counter.Lost = 0;
% for j = 1: size(Content,1)
%
%     if mod(j,1000) == 0
%         clc; disp(['Condition study   - row:', num2str(j)]);
%     end
%
%     %chip top:
%     if isempty(Content{j, column.ChipTop})
%         counter.ChipTop = counter.ChipTop +1;
%     end
%
%     %chip bottom:
%     if isempty(Content{j, column.ChipBot})
%         counter.ChipBot = counter.ChipBot +1;
%     end
%
%     %cont top:
%     if isempty(Content{j, column.ContTop})
%         counter.ContTop = counter.ContTop +1;
%     end
%
%     %cont bottom:
%     if isempty(Content{j, column.ContBot})
%         counter.ContBot = counter.ContBot +1;
%     end
%
%     %Lost:
% %     if ~isempty(Content{j, column.Lost}) && strcmp(Content{j, column.Lost}, 'Lost')
% %         counter.ContLost = counter.Lost +1;
% %     end
% end

rate.ChipTop = counter.ChipTop / (size(Text,1)-1) *100;
rate.ChipBot = counter.ChipBot / (size(Text,1)-1) *100;
rate.ContTop = counter.ContTop / (size(Text,1)-1) *100;
rate.ContBot = counter.ContBot / (size(Text,1)-1) *100;
% rate.Lost = counter.Lost / size(Content,1) *100;

%% display result:

pareto([rate.ChipTop, rate.ChipBot, rate.ContTop, rate.ContBot]);

%% display the run time:
disp(['Elapsed time:', num2str(toc)]);
