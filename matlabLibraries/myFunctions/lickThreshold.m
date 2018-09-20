function timeStampsLick = lickThreshold(animalName, sessionName)

% [timeStampsLick_raw, samplesLick] = Nlx2MatCSC(strcat('D:\Locus Coeruleus Project\Data\', animalName, '\', sessionName, ...
%         '\Unsorted\CSC33.ncs'), [1 0 0 0 1], 0, 1, 1);
[timeStampsLick_raw, samplesLick] = Nlx2MatCSC(strcat('Z:\', animalName, '\', sessionName, ...
        '\Unsorted\CSC33.ncs'), [1 0 0 0 1], 0, 1, 1);
% diffSamplesLick = abs(min(diff(samplesLick))); 
% threshold = 20;
% minimumLickTime = 65; % minimum time between licks (to avoid over counting)
% temp_indices = diffSamplesLick >= threshold;


%% Use the following two lines when lick data isn't perfect
diffSamplesLick = diff(max(samplesLick)); 
threshold = -250; % empirically determined threshold; requires some tweaking to get just right

minimumLickTime = 65; % minimum time between licks (to avoid over counting)
temp_indices = diffSamplesLick <= threshold;

%%

diffSamplesLick = 0;
diffSamplesLick(temp_indices) = 1;
diffSamplesLick = conv(diffSamplesLick, [1 -1]);
diffSamplesLick = diffSamplesLick == 1;
timeStampsLick = timeStampsLick_raw(diffSamplesLick);
timeStampsLick = round(timeStampsLick/1000);
% timeStampsLick(diff(timeStampsLick) < 65) = []; % remove licks that occur too quickly within one another

removeInds = [0 diff([0 diff(timeStampsLick)<minimumLickTime])] == 1;
while any(removeInds)
    timeStampsLick(removeInds) = [];
    removeInds = [0 diff([0 diff(timeStampsLick)<minimumLickTime])] == 1;
end

% plot(timeStampsLick_raw,max(samplesLick)); hold on; plot(timeStampsLick*1000,ones(1,length(timeStampsLick)),'ro');% xlim([6.826 6.836]*10^9); ylim([-1.5 1]*10^4)

% timeStampsLick([0 diff([0 diff(timeStampsLick)<minimumLickTime])] == 1) =
% []; %only remove the first lick that occurs too quickly

% plot(timeStampsLick_raw, mean(samplesLick))
% hold on;
% plot(timeStampsLick*1000,ones(1,length(timeStampsLick)),'ro')