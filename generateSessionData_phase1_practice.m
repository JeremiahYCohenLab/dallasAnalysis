function [sessionData, blockSwitch, blockProbs] = generateSessionData_phase1(~)

% Determine if computer is PC or Mac and set roots and separators appropriately
% [root, sep] = currComputer();
% [animalName, date] = strtok(sessionName, 'd'); 
% animalName = animalName(2:end);
% date = date(1:9);
% sessionFolder = ['m' animalName date];


behavioralDataPath = ['C:\Users\dalla\Desktop\mCG59d20180907.asc'];
sessionText = importData_operantMatching(behavioralDataPath);

sessionData.trialType = [];
sessionData.trialEnd = [];
sessionData.CSon = [];
sessionData.licksL = [];
sessionData.licksR = [];
sessionData.rewardL = [];
sessionData.rewardR = [];
sessionData.rewardTime = [];
sessionData.rewardProbL = [];
sessionData.rewardProbR = [];

blockSwitch = 1;
blockProbs = {};

i = 1;
while isempty(regexp(sessionText{i},'Block Switch at Trial ', 'once'))
    i = i + 1;
end
blockProbs = [blockProbs {sessionText{i}(end-4:end)}];


for i = 1:length(sessionText)
    % determine beginning and end of trial
    if regexp(sessionText{i},'Trial ') == 1 % trial begin 
        temp1 = regexp(sessionText{i},'('); temp2 = regexp(sessionText{i},')');
        currTrial = str2double(sessionText{i}(temp1(1)+1:temp2(1)-1)); % current trial is in between parentheses
        
        tBegin = i; % first index of trial is where the text says 'Trial '
        
        tEndFlag = false;        
        j = i + 1; % start looking for last index of trial
        while (~tEndFlag) 
            if regexp(sessionText{j},'Trial ') == 1
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
            if ~isempty(strfind(sessionText{currTrialInd}, 'Contingency'))
                temp = regexp(sessionText{currTrialInd}, '). ', 'split');
                temp2 = regexp(temp(1,2), '/', 'split');
                sessionData(currTrial).rewardProbL = str2double(temp2{1}{1});
                temp3 = regexp(temp2{1}{2}, ':', 'split');
                sessionData(currTrial).rewardProbR = str2double(temp3(1,1));
            end
            if regexp(sessionText{currTrialInd},'L: ') == 1
                temp = regexp(sessionText(currTrialInd), ': ', 'split');
                allL_licks = [allL_licks str2double(temp{1}{2})];
            elseif regexp(sessionText{currTrialInd},'R: ') == 1
                temp = regexp(sessionText(currTrialInd), ': ', 'split');
                allR_licks = [allR_licks str2double(temp{1}{2})];
            end
            if (~waterDeliverFlag)
                if strfind(sessionText{currTrialInd},'WATER L DELIVERED')
                    temp = regexp(sessionText(currTrialInd), ': ', 'split');
                    sessionData(currTrial).rewardL = 1;
                    sessionData(currTrial).rewardR = NaN;
                    sessionData(currTrial).rewardTime = str2double(temp{1}{2});
                    waterDeliverFlag = true;
                elseif strfind(sessionText{currTrialInd},'WATER L NOT DELIVERED')
                    temp = regexp(sessionText(currTrialInd), ': ', 'split');
                    sessionData(currTrial).rewardL = 0;
                    sessionData(currTrial).rewardR = NaN;
                    sessionData(currTrial).rewardTime = str2double(temp{1}{2});
                    waterDeliverFlag = true;
                elseif strfind(sessionText{currTrialInd},'WATER R DELIVERED')
                    temp = regexp(sessionText(currTrialInd), ': ', 'split');
                    sessionData(currTrial).rewardR = 1;
                    sessionData(currTrial).rewardL = NaN;
                    sessionData(currTrial).rewardTime = str2double(temp{1}{2});
                    waterDeliverFlag = true;
                elseif strfind(sessionText{currTrialInd},'WATER R NOT DELIVERED')
                    temp = regexp(sessionText(currTrialInd), ': ', 'split');
                    sessionData(currTrial).rewardR = 0;
                    sessionData(currTrial).rewardL = NaN;
                    sessionData(currTrial).rewardTime = str2double(temp{1}{2});
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
                sessionData(currTrial).licksL = allL_licks;
                sessionData(currTrial).licksR = allR_licks;
                if ~waterDeliverFlag
                    sessionData(currTrial).rewardL = NaN;
                    sessionData(currTrial).rewardR = NaN;
                    sessionData(currTrial).rewardTime = NaN;
                end
                if tEnd ~= length(sessionText)
                    temp = regexp(sessionText(tEnd+1), ': ', 'split');
                    sessionData(currTrial).trialEnd = str2double(temp{1}{2});
                else
                    sessionData(currTrial).trialEnd = NaN;
                end
            end
        end
        if regexp(sessionText{currTrialInd},'Block Switch at Trial ') == 1
            if currTrial ~= 1
                blockSwitch = [blockSwitch currTrial];
                blockProbs = [blockProbs {sessionText{currTrialInd}(end-4:end)}];
            end
        end
    end
end



% savepath = [behavioralDataPath(1:strfind(behavioralDataPath,'behavior')-1) 'sorted' sep 'session' sep];
% if isempty(dir(savepath))
%     mkdir(savepath)
end

% f_IndA = find(behavioralDataPath==sep,1,'last');
% f_IndB = strfind(behavioralDataPath,'.asc');
% filename = behavioralDataPath(f_IndA+1:f_IndB-1);

% save([savepath filename '_sessionData_behav.mat'], 'sessionData', 'blockSwitch', 'blockProbs');
% end

function dataOutput = importData_operantMatching(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   MBB039D20160712 = IMPORTFILE(FILENAME) Reads data from text file
%   FILENAME for the default selection.
%
%   MBB039D20160712 = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data
%   from rows STARTROW through ENDROW of text file FILENAME.
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
fileID = fopen(filename,'r');

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