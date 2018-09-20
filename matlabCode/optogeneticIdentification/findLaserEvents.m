function [laserStruct] = findLaserEvents(filename)

    %% Import Event TTLs, Event TTL Timestamps, and Cluster Data
    % Import Event TTLs and Timestamps
    
    [root, sep] = currComputer();

    [animalName, date] = strtok(filename, 'd'); 
    animalName = animalName(2:end);
    date = date(1:9);
    sessionFolder = ['m' animalName date];
    sessionDataPath = [root animalName sep sessionFolder sep 'ephys' sep 'opto' sep 'Events.nev'];

    [timeStampsTTLs, TTLs] = Nlx2MatEV(sessionDataPath, [1 0 1 0 0], 0, 1);
    events = [timeStampsTTLs; TTLs];
    
     % Import cluster data
    sortedSpikesPath = [root animalName sep sessionFolder sep 'sorted' sep 'opto' sep];
    allFiles = dir(strcat(sortedSpikesPath));
    for j = find(~cellfun(@isempty,strfind({allFiles.name},'TT')) & ~cellfun(@isempty,strfind({allFiles.name},'txt'))); % load all clustered spikes
        load(strcat(sortedSpikesPath, allFiles(j).name))
    end

    %% Parse Out Laser Time Data
    offset = 0; % all offset
    laser = 64; % laser on
    
%     If laser is present, find all laser on times and laser off times
    laserStruct = [];
    eventsLaser = [];
    errorInds = [];
    if any(events(2,:) >= laser)
        errorCount = 0;
        for i = 1:length(events)
            if ~errorCount
                if events(2,i) >= laser                         %gets rid of sequences where laser is interrupted by other signal
                    foundFlag = 0;                      
                    tmpInd = 1;
                    while ~foundFlag
                        if events(2,i) >= laser & events(2,i+tmpInd) >= laser
                            tmpInd = tmpInd+1;
                            errorCount = errorCount+1;
                        else
                            eventsLaser = [eventsLaser events(:,i) events(:,i+tmpInd)];
                            foundFlag = 1;
                        end
                    end
                end
            else
                errorCount = errorCount -1;
            end
        end
        for i = 1:length(eventsLaser)-1                             %finds all errors due to conflicting signals
            if eventsLaser(1,i) == eventsLaser(1,i+1)
                errorInds = [errorInds i];
            end
        end
        if ~isempty(errorInds)                                      %removes erroneous signals
            for i=1:length(errorInds)
                eventsLaser = [eventsLaser(:,1:(errorInds(i)-1)) eventsLaser(:,(errorInds(i)+2):end)];
            end
        end
        
        temp1 = diff(eventsLaser(1,:)); % all time differences
        laserOn = find(eventsLaser(2,:) >= laser); % all laser on indices
        laserHist = temp1(laserOn); % all laser on times differences
        laserStruct.pulseWidth = mode(laserHist); % mode of above distribution is the laser time
        
        
        laserStruct.laserOn = eventsLaser(1,laserOn(laserHist >= laserStruct.pulseWidth-100 & laserHist < laserStruct.pulseWidth+100)); % all valid laser on times based on pulse width
        laserStruct.laserOff = eventsLaser(1,1+laserOn(laserHist >= laserStruct.pulseWidth-100 & laserHist < laserStruct.pulseWidth+100)); % all valid laser off times based on pulse width
        
        troubleInds = find(~(laserHist >= laserStruct.pulseWidth-100 & laserHist < laserStruct.pulseWidth+100)); % indices with an inappropriate laser off time
 
        for j = 1:length(troubleInds)
            tempInd = 2;
            foundFlag = 0;
            while ~foundFlag
                trouble_laserOn = eventsLaser(1,laserOn(troubleInds(j)));
                trouble_laserOff = eventsLaser(1,tempInd+laserOn(troubleInds(j)));
                foundFlag = trouble_laserOff - trouble_laserOn >= laserStruct.pulseWidth - 100 & trouble_laserOff - trouble_laserOn < laserStruct.pulseWidth + 100; % if the appropriate index is found
                if trouble_laserOff - trouble_laserOn >= laserStruct.pulseWidth + 100
                    trouble_laserOn = [];
                    trouble_laserOff = [];
                    foundFlag = 1;
                end
                tempInd = tempInd + 1;
            end
            laserStruct.laserOn = [laserStruct.laserOn trouble_laserOn];
            laserStruct.laserOff = [laserStruct.laserOff trouble_laserOff];
        end
        laserStruct.laserOn = sort(laserStruct.laserOn/1000); %sort and convert from us to ms
        laserStruct.laserOff = sort(laserStruct.laserOff/1000); %sort and convert from us to ms
        
    end
    
        
        allVar = who;
        allClust = find(~cellfun(@isempty,strfind(allVar,'SS')) & ~cellfun(@isempty,strfind(allVar,'TT')))';
        for m = allClust
            eval(sprintf('laserStruct.%s = %s''/1000;', allVar{m}, allVar{m}))
        end
        
    save(strcat(sortedSpikesPath, 'laserStruct.mat'), 'laserStruct')
end
