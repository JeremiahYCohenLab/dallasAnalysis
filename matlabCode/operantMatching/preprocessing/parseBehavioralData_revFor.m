function outputStruct = parseBehavioralData_revFor(sessionData, blockSwitch)

responseInds = find(~isnan([sessionData.rewardTime]));
omitInds = isnan([sessionData.rewardTime]);

tempBlockSwitch = blockSwitch;
for i = 2:length(blockSwitch)
    subVal = sum(omitInds(tempBlockSwitch(i-1):tempBlockSwitch(i)));
    blockSwitch(i:end) = blockSwitch(i:end) - subVal;
end

allReward_R = [sessionData(responseInds).rewardR]; 
allReward_L = [sessionData(responseInds).rewardL]; 
allChoices = NaN(1,length(sessionData(responseInds)));
allChoices(~isnan(allReward_R)) = 1;
allChoices(~isnan(allReward_L)) = -1;

allReward_R(isnan(allReward_R)) = 0;
allReward_L(isnan(allReward_L)) = 0;
allChoice_R = double(allChoices == 1);
allChoice_L = double(allChoices == -1);

allRewards = zeros(1,length(allChoices));
allRewards(logical(allReward_R)) = 1;
allRewards(logical(allReward_L)) = -1;

timeBtwn = [[sessionData(2:end).rewardTime] - [sessionData(1:end-1).CSon]];
timeBtwn(timeBtwn < 0 ) = 0;
timeBtwn = [0 timeBtwn(~isnan(timeBtwn))/1000];

%% Convert to outputStruct
outputStruct.responseInds = responseInds;
outputStruct.omitInds = omitInds;
outputStruct.blockSwitch = blockSwitch;
outputStruct.allRewards = allRewards;
outputStruct.allReward_R = allReward_R;
outputStruct.allReward_L = allReward_L;
outputStruct.allChoices = allChoices;
outputStruct.allChoice_R = allChoice_R;
outputStruct.allChoice_L = allChoice_L;
outputStruct.timeBtwn = timeBtwn;