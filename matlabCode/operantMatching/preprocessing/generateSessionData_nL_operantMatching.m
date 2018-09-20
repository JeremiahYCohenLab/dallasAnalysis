function [sessionData] = generateSessionData_nL_operantMatching(sessionName)
%generateSessionData_neuraLynx_operantMatching will convert Neuralynx data to a structure
%   [sessionData] = generateSessionData_neuraLynx_operantMatching(sessionFolder)
%   INPUT: sessionFolder
%       e.g. sessionFolder = 'mBB040d20160930'
%       Assumes there is a 'neuralynx' and 'sorted' folder within sessionFolder
%       Neuralynx data must contain behavioral data (in the form of TTL inputs in Events.nev)
%   OUTPUT: sessionData


% Determine if computer is PC or Mac and set roots and separators appropriately
[root,sep] = currComputer();

%% Assign variables
% Path
[animalName] = strtok(sessionName, 'd');
animalName = animalName(2:end);
if isstrprop(sessionName(end), 'alpha')
    filepath = [root animalName sep sessionName(1:end-1) sep 'ephys' sep 'session ' sessionName(end) sep];
    sortedDir = [root animalName sep sessionName(1:end-1) sep 'sorted' sep 'session ' sessionName(end) sep];
else
    filepath = [root animalName sep sessionName sep 'ephys' sep 'session' sep];
    sortedDir = [root animalName sep sessionName sep 'sorted' sep 'session' sep];
end

% TTL values
CSplus = 136;
CSminus = 160;
waterR = 1; % pin 10
waterL = 2; % pin 11
lickR = 32; % pin 15
lickL = 16; % pin 14
offset = 0;

% Notes
%   odorBlank = 8 if preceding 136;
%   odorBlank = 128 if immediately turning off (precedes a 0)
%   lickR + waterR = 33;
%   lickL + waterL = 18;

%% Import event TTLs, TTL timestamps, and cluster data (if it exists)
% TTLs
[timestamps, TTL] = Nlx2MatEV([filepath 'Events.nev'],[1 0 1 0 0], 0, 1);
ttlEvents = [timestamps; TTL];

% Clustered data
sortedFiles = dir(sortedDir);
for j = find(~cellfun(@isempty,strfind({sortedFiles.name},'TT')) & ...
             ~cellfun(@isempty,strfind({sortedFiles.name},'txt')))
    load([sortedDir sortedFiles(j).name]);
end

%% Generate structure
ttlEvents(1,:) = round(ttlEvents(1,:)/1000); % convert from us to ms
ttlEvents = ttlEvents(:,ttlEvents(2,:) == CSplus | ttlEvents(2,:) ==  CSminus | ttlEvents(2,:) == waterR | ...
                        ttlEvents(2,:) == waterL | ttlEvents(2,:) == lickR | ttlEvents(2,:) == lickL | ...
                        ttlEvents(2,:) == (lickR+waterR) | ttlEvents(2,:) == (lickL+waterL) | ...
                        ttlEvents(2,:) == offset); % remove garbage data
% allOffs = ttlEvents(2,:) == offset; allOffs = allOffs(:,1:end-1);
% ttlEvents(:, find(allOffs & diff(ttlEvents(2,:)) == 0) + 1) = [];

% Initialize structure
sessionData.trialType = [];
sessionData.trialEnd = [];
sessionData.CSon = [];
sessionData.licksL = [];
sessionData.licksR = [];
sessionData.rewardL = [];
sessionData.rewardR = [];
sessionData.rewardTime = [];
sessionData.allSpikes = [];

lickWindowDuration = 2000;

% Find tetrode data (files with SS and TT in the name)
allVar = who;
allClust = find(~cellfun(@isempty,strfind(allVar,'SS')) & ~cellfun(@isempty,strfind(allVar,'TT')))';
for m = allClust
    eval(sprintf('%s = round(%s/1000);', allVar{m}, allVar{m})) % convert from us to ms
    sessionData(m-2).allSpikes = transpose(eval(allVar{m}));
end


trial = 1;
for i = 1:length(ttlEvents)
    if ttlEvents(2,i) == CSplus || ttlEvents(2,i) == CSminus % if this is a trial begin and continue to end of trial
        if ttlEvents(2,i) == CSplus
            sessionData(trial).trialType = 'CSplus';
        elseif ttlEvents(2,i) == CSminus
            sessionData(trial).trialType = 'CSminus';
        end
        sessionData(trial).CSon = ttlEvents(1,i);
        sessionData(trial).rewardR = NaN; 
        sessionData(trial).rewardL = NaN;
        sessionData(trial).rewardTime = NaN;
        
        currInd = i+1; % begin at the next index
        tEndFlag = false; % end when this is true
        while tEndFlag == false
            if ttlEvents(2,currInd) == lickR % right licks
                sessionData(trial).licksR = [sessionData(trial).licksR ttlEvents(1,currInd)];
            elseif ttlEvents(2,currInd) == lickL % left licks
                sessionData(trial).licksL = [sessionData(trial).licksL ttlEvents(1,currInd)];
            elseif ttlEvents(2,currInd) == (lickR+waterR) % R rewards; waterR by itself could be evap code
                if isnan(sessionData(trial).rewardR)
                    sessionData(trial).rewardR = 1;
                    sessionData(trial).rewardTime = ttlEvents(1,currInd);
                end
            elseif ttlEvents(2,currInd) == (lickL+waterL) % L rewards; waterL by itself could be evap code
                if isnan(sessionData(trial).rewardL)
                    sessionData(trial).rewardL = 1;
                    sessionData(trial).rewardTime = ttlEvents(1,currInd);
                end
            end
            
            currInd = currInd + 1; % continue to the next index
            if currInd <= length(ttlEvents) % if not at the end of the file, stop when index is on the next trial
                if ttlEvents(2,currInd) == CSplus || ttlEvents(2,currInd) == CSminus
                    sessionData(trial).trialEnd = ttlEvents(1,currInd);
                    tEndFlag = true;
                end
            else % if at the end of the file, escape
                sessionData(trial).trialEnd = NaN;
                tEndFlag = true;
            end
        end
        
        % whole bunch of code to figure out when to indicate no reward given
        if strcmp(sessionData(trial).trialType,'CSplus') % only makes sense on CS+ trials
            if isnan(sessionData(trial).rewardR) && isnan(sessionData(trial).rewardL) % if fields are default
                if ~isempty(sessionData(trial).licksR) && ~isempty(sessionData(trial).licksL) % L and R licks
                    lickRlat = sessionData(trial).licksR(1) - sessionData(trial).CSon;
                    lickLlat = sessionData(trial).licksL(1) - sessionData(trial).CSon;
                    if lickRlat <= lickWindowDuration && lickLlat <= lickWindowDuration % L and R within lick window; shouldn't ever happen
                        if lickRlat < lickLlat
                            sessionData(trial).rewardR = 0;
                            sessionData(trial).rewardTime = sessionData(trial).licksR(1);
                        else
                            sessionData(trial).rewardL = 0;
                            sessionData(trial).rewardTime = sessionData(trial).licksL(1);
                        end
                    elseif lickRlat <= lickWindowDuration % R lick within lick window
                        sessionData(trial).rewardR = 0;
                        sessionData(trial).rewardTime = sessionData(trial).licksR(1);
                    elseif lickLlat <= lickWindowDuration % L lick within lick window
                        sessionData(trial).rewardL = 0;
                        sessionData(trial).rewardTime = sessionData(trial).licksL(1);
                    end
                elseif ~isempty(sessionData(trial).licksR) && isempty(sessionData(trial).licksL) % R lick window
                    lickRlat = sessionData(trial).licksR(1) - sessionData(trial).CSon;
                    if lickRlat <= lickWindowDuration % R lick within lick window
                        sessionData(trial).rewardR = 0;
                        sessionData(trial).rewardTime = sessionData(trial).licksR(1);
                    end
                elseif isempty(sessionData(trial).licksR) && ~isempty(sessionData(trial).licksL) % L lick window
                    lickLlat = sessionData(trial).licksL(1) - sessionData(trial).CSon;
                    if lickLlat <= lickWindowDuration % L lick within lick window
                        sessionData(trial).rewardL = 0;
                        sessionData(trial).rewardTime = sessionData(trial).licksL(1);
                    end
                end
            end
        end
        
        % remove bounced lick signals (those that occur within 55ms of one another)
        if ~isempty(sessionData(trial).licksR)
            sessionData(trial).licksR = removeBadLicks(sessionData(trial).licksR);
        end
        if ~isempty(sessionData(trial).licksL)
            sessionData(trial).licksL = removeBadLicks(sessionData(trial).licksL);
        end
        
        % add spike data
        L = sessionData(trial).CSon;
        U = sessionData(trial).trialEnd;
        if isnan(U)
            U = ttlEvents(1,end);
        end
        for l = allClust
            eval(sprintf('sessionData(trial).%s = %s(%s >= L & %s < U)'';', ...
            allVar{l}, allVar{l}, allVar{l}, allVar{l}))
        end
        
        trial = trial+1;
    end
end

% save the data
if isempty(dir(sortedDir))
    mkdir(sortedDir)
end

% sessionFolder(~(sessionFolder==sep)) removes the separator (/ or \) in the filename before saving
save([sortedDir sep sessionName(~(sessionName==sep)) '_sessionData_nL.mat'], 'sessionData')