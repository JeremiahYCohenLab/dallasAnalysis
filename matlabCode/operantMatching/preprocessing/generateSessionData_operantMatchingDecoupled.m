function [behSessionData, blockSwitch, blockSwitchL, blockSwitchR] = generateSessionData_operantMatchingDecoupled(sessionName)


% Determine if computer is PC or Mac and set roots and separators appropriately
[root, sep] = currComputer();

% Generate the correct file path and see what files are available
[animalName, date] = strtok(sessionName, 'd'); 
animalName = animalName(2:end);
date = date(1:9);
sessionFolder = ['m' animalName date];
filepath = [root animalName sep sessionFolder sep 'behavior' sep sessionName '.asc'];

if isstrprop(sessionName(end), 'alpha')
    savepath = [root animalName sep sessionFolder sep 'sorted' sep 'session ' sessionName(end) sep];
else
    savepath = [root animalName sep sessionFolder sep 'sorted' sep 'session' sep];
end

sessionText = importData_operantMatching(filepath);

behSessionData.trialType = [];
behSessionData.trialEnd = [];
behSessionData.CSon = [];
behSessionData.licksL = [];
behSessionData.licksR = [];
behSessionData.rewardL = [];
behSessionData.rewardR = [];
behSessionData.rewardTime = [];
behSessionData.rewardProbL = [];
behSessionData.rewardProbR = [];
behSessionData.allLicks = [];

blockSwitch = 1;
blockSwitchL = 1;
blockSwitchR = 1;

for i = 1:length(sessionText)
    % determine beginning and end of trial
    if regexp(sessionText{i},'L Trial ') == 1 % trial begin 
        temp1 = regexp(sessionText{i},'('); temp2 = regexp(sessionText{i},')');
        currTrial = str2double(sessionText{i}(temp1(1)+1:temp2(1)-1)); % current trial is in between parentheses
        
        tBegin = i; % first index of trial is where the text says 'Trial '
        
        tEndFlag = false;        
        j = i + 1; % start looking for last index of trial
        while (~tEndFlag) 
            if regexp(sessionText{j},'L Trial ') == 1
                tEnd = j - 1; 
                tEndFlag = true;
            else
                j = j + 1;
                if j == length(sessionText)
                    tEnd = length(sessionText);
                    tEndFlag = true;
                end
            end
        end
        
        waterDeliverFlag = false;
        allL_licks = [];
        allR_licks = [];
        for currTrialInd = tBegin+1:tEnd
            if strfind(sessionText{currTrialInd}, 'Contingency') == 1
                temp = regexp(sessionText{currTrialInd}, '). ', 'split');
                temp2 = regexp(temp(1,2), '/', 'split');
                behSessionData(currTrial).rewardProbL = str2double(temp2{1}{1});
                temp3 = regexp(temp2{1}{2}, ':', 'split');
                behSessionData(currTrial).rewardProbR = str2double(temp3(1,1));
            end
            if strfind(sessionText{currTrialInd},'CS PLUS')
                behSessionData(currTrial).trialType = 'CSplus';
                temp = regexp(sessionText(currTrialInd), ': ', 'split');
                behSessionData(currTrial).CSon = str2double(temp{1}{2});
            elseif strfind(sessionText{currTrialInd},'CS MINUS')
                behSessionData(currTrial).trialType = 'CSminus';
                temp = regexp(sessionText(currTrialInd), ': ', 'split');
                behSessionData(currTrial).CSon = str2double(temp{1}{2});
            end
            if regexp(sessionText{currTrialInd},'L: ') == 1
                temp = regexp(sessionText(currTrialInd), ': ', 'split');
                allL_licks = [allL_licks str2double(temp{1}{2})];
                behSessionData(1).allLicks = [behSessionData(1).allLicks str2double(temp{1}{2})];
            elseif regexp(sessionText{currTrialInd},'R: ') == 1
                temp = regexp(sessionText(currTrialInd), ': ', 'split');
                allR_licks = [allR_licks str2double(temp{1}{2})];
                behSessionData(1).allLicks = [behSessionData(1).allLicks str2double(temp{1}{2})];
            end
            if (~waterDeliverFlag)
                if strfind(sessionText{currTrialInd},'WATER L DELIVERED')
                    temp = regexp(sessionText(currTrialInd), ': ', 'split');
                    behSessionData(currTrial).rewardL = 1;
                    behSessionData(currTrial).rewardR = NaN;
                    behSessionData(currTrial).rewardTime = str2double(temp{1}{2});
                    waterDeliverFlag = true;
                elseif strfind(sessionText{currTrialInd},'WATER L NOT DELIVERED')
                    temp = regexp(sessionText(currTrialInd), ': ', 'split');
                    behSessionData(currTrial).rewardL = 0;
                    behSessionData(currTrial).rewardR = NaN;
                    behSessionData(currTrial).rewardTime = str2double(temp{1}{2});
                    waterDeliverFlag = true;
                elseif strfind(sessionText{currTrialInd},'WATER R DELIVERED')
                    temp = regexp(sessionText(currTrialInd), ': ', 'split');
                    behSessionData(currTrial).rewardR = 1;
                    behSessionData(currTrial).rewardL = NaN;
                    behSessionData(currTrial).rewardTime = str2double(temp{1}{2});
                    waterDeliverFlag = true;
                elseif strfind(sessionText{currTrialInd},'WATER R NOT DELIVERED')
                    temp = regexp(sessionText(currTrialInd), ': ', 'split');
                    behSessionData(currTrial).rewardR = 0;
                    behSessionData(currTrial).rewardL = NaN;
                    behSessionData(currTrial).rewardTime = str2double(temp{1}{2});
                    waterDeliverFlag = true;
                end
            end
%             if strfind(sessionText{currTrialInd},'AUTOMATIC WATER L DELIVERED')
%                 temp = regexp(sessionText(currTrialInd), ': ', 'split');
%                 sessionData(currTrial).autoWaterL = str2double(temp{1}{2});
%             elseif strfind(sessionText{currTrialInd},'AUTOMATIC WATER R DELIVERED')
%                 temp = regexp(sessionText(currTrialInd), ': ', 'split');
%                 sessionData(currTrial).autoWaterR = str2double(temp{1}{2});
%             end
            
            if currTrialInd == tEnd % run this at the last index || currTrialInd == length(sessionText)-1
                behSessionData(currTrial).licksL = allL_licks;
                behSessionData(currTrial).licksR = allR_licks;
                if ~waterDeliverFlag
                    behSessionData(currTrial).rewardL = NaN;
                    behSessionData(currTrial).rewardR = NaN;
                    behSessionData(currTrial).rewardTime = NaN;
                end
                if tEnd ~= length(sessionText)
                    temp = regexp(sessionText(tEnd+3), ': ', 'split');
                    behSessionData(currTrial).trialEnd = str2double(temp{1}{2});
                else
                    behSessionData(currTrial).trialEnd = NaN;
                end
            end
        end
        if regexp(sessionText{currTrialInd},'L Block Switch at Trial ') == 1
            if currTrial ~= 1
                blockSwitch = [blockSwitch currTrial];
                blockSwitchL = [blockSwitchL currTrial];
            end
        end
        if regexp(sessionText{currTrialInd},'R Block Switch at Trial ') == 1
            if currTrial ~= 1
                blockSwitch = [blockSwitch currTrial];
                blockSwitchR = [blockSwitchR currTrial];
            end
        end
    end
end


if isempty(dir(savepath))
    mkdir(savepath)
end
save([savepath sessionName '_sessionData_behav.mat'], 'behSessionData', 'blockSwitch', 'blockSwitchL', 'blockSwitchR');
end

function dataOutput = importData_operantMatching(sessionName, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   MBB039D20160712 = IMPORTFILE(sessionName) Reads data from text file
%   sessionName for the default selection.
%
%   MBB039D20160712 = IMPORTFILE(sessionName, STARTROW, ENDROW) Reads data
%   from rows STARTROW through ENDROW of text file sessionName.
%
% Example:
%   mBB039d20160712 = importfile('mBB039d20160712.asc', 1, 5986);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2016/07/12 16:31:53

%% Initialize variables.
delimiter = '';
if nargin<=2
    startRow = 1;
    endRow = inf;
end

%% Format string for each line of text:
%   column1: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%[^\n\r]';

%% Open the text file.
fileID = fopen(sessionName,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'ReturnOnError', false);
    dataArray{1} = [dataArray{1};dataArrayBlock{1}];
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
dataOutput = [dataArray{1:end-1}];
end