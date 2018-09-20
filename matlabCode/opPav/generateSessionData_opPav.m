function [sessionData, blockSwitch] = generateSessionData_opPav(filename)

% Determine if computer is PC or Mac and set roots and separators appropriately
[root, sep] = currComputer();

% Generate the correct file path and see what files are available
[animalName, date] = strtok(filename, 'd'); 
animalName = animalName(2:end);
date = date(1:9);
sessionFolder = ['m' animalName date];
filepath = [root animalName sep sessionFolder sep 'behavior' sep filename '.asc'];

sessionText = importData_opPav(filepath);

sessionData.trialType = [];
sessionData.blockType = [];
sessionData.trialIndex = [];
sessionData.trialEnd = [];
sessionData.CSon = [];
sessionData.lick = [];
sessionData.pull = [];
sessionData.reward = [];
sessionData.rewardTime = [];
sessionData.airpuff = [];
sessionData.airpuffTime = [];
sessionData.avoidance = [];
sessionData.avoidanceTime = [];

blockData.blockNumber = [];
blockData.blockType = [];
blockData.blockStart = [];
blockData.blockSwitch = [];

allLicks = [];
allPulls = [];
allManualWater = [];
currBlock = 0;
blockSwitchFlag = false;


for i = 1:length(sessionText)
    
    if strfind(sessionText{i}, 'Block Number')          %find block number
        currBlock = currBlock +1;
        blockData(currBlock).blockNumber = currBlock;
    end
           
    if strfind(sessionText{i},'Block Start')            %block start time
       temp = regexp(sessionText(i), ': ', 'split');
       blockData(currBlock).blockStart = str2double(temp{1}{2});
       blockSwitchFlag = true;
    end
    if strfind(sessionText{i},'Block Type')              %block type
        temp = regexp(sessionText(i), ': ', 'split');
        blockData(currBlock).blockType = temp{1}{2};
    end
    if regexp(sessionText{i},'lick') == 1                %find all licks
        temp = regexp(sessionText(i), ': ', 'split');
        allLicks = [allLicks str2double(temp{1}{2})];
    end
    if regexp(sessionText{i},'joystick') == 1            %find all pulls
        temp = regexp(sessionText(i), ': ', 'split');
        allPulls = [allPulls str2double(temp{1}{2})];
    end
    if regexp(sessionText{i},'manual') == 1              %find all manual water deliveries
        temp = regexp(sessionText(i), ': ', 'split');
        allManualWater = [allPulls str2double(temp{1}{2})];
    end

    
    % determine beginning and end of trial
    if regexp(sessionText{i},'Trial:') == 1 % trial begin
        if strfind(sessionText{i}, '|')
            temp1 = regexp(sessionText{i},': '); temp2 = regexp(sessionText{i},'     |');
            currTrial = str2double(sessionText{i}(temp1+2:temp2-1)); % current trial is in between parentheses
        else
            temp = regexp(sessionText{i}, ':', 'split');        % MIGHT NEED TO CHANGE W NEW BEH CODE
            currTrial = temp{2};
        end
        
        if blockSwitchFlag == true
            blockData(currBlock).blockSwitch = currTrial;
            blockSwitchFlag = false;
        end
            
        
        tBegin = i - 1; % first index of trial is where the text says 'C S +/- '
        
        tEndFlag = false;        
        j = i + 1; % start looking for last index of trial
        while (~tEndFlag) 
            if regexp(sessionText{j},'ITI') == 1
                tEnd = j; 
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
        airpuffDeliverFlag = false;
        airpuffAvoidFlag = false;
        USdeliverFlag = false;
        allLicksTemp = [];
        allPullsTemp = [];
        
        for currTrialInd = tBegin:tEnd
            if strfind(sessionText{currTrialInd},'C S +')                       %find CS type
                sessionData(currTrial).trialType = 'CSplus';
            elseif strfind(sessionText{currTrialInd},'C S -')
                sessionData(currTrial).trialType = 'CSminus';
            end
            
            if strfind(sessionText{currTrialInd}, 'Odor On')
                temp = regexp(sessionText(currTrialInd), ': ', 'split');
                sessionData(currTrial).CSon = str2double(temp{1}{2});
            end
            
            if strfind(sessionText{currTrialInd},'operant reward')          %find block type
                sessionData(currTrial).blockType = 'opRew';
            elseif strfind(sessionText{currTrialInd},'operant airpuff')
                sessionData(currTrial).blockType = 'opAve';
            elseif strfind(sessionText{currTrialInd},'pavlovian reward')
                sessionData(currTrial).blockType = 'pavRew';
            elseif strfind(sessionText{currTrialInd},'pavlovian airpuff')
                sessionData(currTrial).blockType = 'pavAve';
            else
                sessionData(currTrial).blockType = NaN;
            end
            
            if regexp(sessionText{currTrialInd},'lick') == 1                %find all licks
                temp = regexp(sessionText(currTrialInd), ': ', 'split');
                allLicks = [allLicks str2double(temp{1}{2})];
                allLicksTemp = [allLicksTemp str2double(temp{1}{2})];
            end
            if regexp(sessionText{currTrialInd},'joystick') == 1            %find all pulls
                temp = regexp(sessionText(currTrialInd), ': ', 'split');
                allPulls = [allPulls str2double(temp{1}{2})];
                allPullsTemp = [allPullsTemp str2double(temp{1}{2})];
            end
            if regexp(sessionText{currTrialInd},'manual') == 1            %find all manual water deliveries
                temp = regexp(sessionText(currTrialInd), ': ', 'split');
                allManualWater = [allPulls str2double(temp{1}{2})];
            end
            
            
            if (~USdeliverFlag) 
                if strfind(sessionText{currTrialInd},'WATER DELIVERED')
                    temp = regexp(sessionText(currTrialInd), ': ', 'split');
                    sessionData(currTrial).reward = 1;
                    sessionData(currTrial).rewardTime = str2double(temp{1}{2});
                    waterDeliverFlag = true;
                    USdeliverFlag = true;
                elseif strfind(sessionText{currTrialInd},'AIRPUFF DELIVERED')
                    temp = regexp(sessionText(currTrialInd), ': ', 'split');
                    sessionData(currTrial).airpuff = 1;
                    sessionData(currTrial).airpuffTime = str2double(temp{1}{2});
                    airpuffDeliverFlag = true;
                    USdeliverFlag = true;
                elseif strfind(sessionText{currTrialInd},'OP AIRPUFF AVOIDED')
                    temp = regexp(sessionText(currTrialInd), ': ', 'split');
                    sessionData(currTrial).avoidance = 1;
                    sessionData(currTrial).avoidanceTime = str2double(temp{1}{2});
                    airpuffAvoidFlag = true;
                    USdeliverFlag = true;
                end
            end

            
            if currTrialInd == tEnd % run this at the last index || currTrialInd == length(sessionText)-1
                sessionData(currTrial).lick = allLicksTemp;
                sessionData(currTrial).pull = allPullsTemp;
                
                if ~waterDeliverFlag
                    sessionData(currTrial).reward = NaN;
                    sessionData(currTrial).rewardTime = NaN;
                end
                if ~airpuffDeliverFlag
                    sessionData(currTrial).airpuff = NaN;
                    sessionData(currTrial).airpuffTime = NaN;
                end
                if ~airpuffAvoidFlag
                    sessionData(currTrial).avoidance = NaN;
                    sessionData(currTrial).avoidanceTime = NaN;
                end
                if ~USdeliverFlag
                    sessionData(currTrial).reward = NaN;
                    sessionData(currTrial).rewardTime = NaN;
                    sessionData(currTrial).airpuff = NaN;
                    sessionData(currTrial).airpuffTime = NaN;
                    sessionData(currTrial).avoidance = NaN;
                    sessionData(currTrial).avoidanceTime = NaN;
                end
                                        
                
                if tEnd ~= length(sessionText)
                    temp = regexp(sessionText(tEnd), '|   ', 'split');
                    sessionData(currTrial).trialEnd = str2double(temp{1}{3});
                else
                    sessionData(currTrial).trialEnd = NaN;
                end
            end
        end
    end
end

savepath = [root animalName sep sessionFolder sep 'sorted' sep];
if isempty(dir(savepath))
    mkdir(savepath)
end
save([savepath filename '_sessionData.mat'], 'sessionData', 'blockData', 'allLicks', 'allPulls', 'allManualWater');
end

function dataOutput = importData_opPav(filename, startRow, endRow)
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