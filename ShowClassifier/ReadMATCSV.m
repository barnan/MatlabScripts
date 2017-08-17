% Read and convert csv file into a MATLAB structure
% [ data, header ] = ReadMATCSV ( fileName )
%
% fileName  - in: full path of the input file
% header    - out: header of the csv file (this is the colun name, the name of the features)
% data      - out: data in tabular from
%
% Barna N 2016
%
%%
% read in the given mat file or csv file and order the data into the same 
% header and data structure
function [data, header] = ReadMATCSV (fileName)

if strcmp(fileName(end-3:end),'.csv')
    data = ReadCSV(fileName);
    
    header = data(2:3,:);
    data(1,:) = [];
    data(1,:) = [];
    data(1,:) = [];
        
elseif strcmp(fileName(end-3:end),'.mat')
    [data, header] = ReadMAT(fileName);
end

end

% read in the given csv file and order the data in tabular form
function data = ReadCSV (fileName)
try
    data = {};
    fidOfFile = fopen(fileName);
    while ~feof(fidOfFile)
        
        %read in a single line:
        line = fgetl(fidOfFile);
        lineArr = strread(line,'%s','delimiter',';');
        
        %if size(lineArr,1) > 10    % magic number!!
            if size(lineArr,1) < size(data,1)
                lineArr{size(data,1),1} = '';
            end
            if size(lineArr,1) > size(data,1)
                for i = 1: size(data,2)
                    data{size(lineArr,1),i} = '';
                end
            end
            data = [data lineArr];
        %end
    end
    data = data';
    fclose(fidOfFile);
catch exc
    msg = ['Problem during CSV file reading: ',fileName,' '];
    warning(msg);
    rethrow (exc);
end
end

% Read in the given mat file
function [data, header] = ReadMAT (fileName)

if ~strcmp(fileName(end-3:end),'.mat')
    fileName = [fileName '.mat'];
else
    %fileName = fileName(1:end-4);
end

try
    position = strfind(fileName,'\');
    variabName = fileName(position(end)+1:end-4);
    
    eval (['load ' fileName]);
    
    variabList = who;
    
    if find(strcmp(variabList,'Header'))
        eval('header = Header;');
    elseif find(strcmp(variabList,'Header_full'))
        eval('header = Header_full;');
    else
        header = {};
    end
        
    if find(strcmp(variabList,variabName))
        eval(['data = ',variabName,';']);
    else
        data = [];
    end
catch exc
    msg = ['Problem during MAT file reading: ',fileName,' '];
    warning(msg);
    rethrow (exc);
end
end



