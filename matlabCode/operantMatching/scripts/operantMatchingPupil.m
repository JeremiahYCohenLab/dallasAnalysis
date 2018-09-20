sessionPupilData = sessionPupilData(1:end-1);
% sessionPupilData(284:336) = [];

allReward_R = [sessionPupilData.rewardR];
allReward_L = [sessionPupilData.rewardL];
allChoice = NaN(1,length(sessionPupilData));
allChoice(~isnan(allReward_R)) = 1;
allChoice(~isnan(allReward_L)) = -1;

relevInds = find(~isnan(allChoice));

% for i = 284:336
% relevInds(find(relevInds == i)) = [];
% end; relevInds(end) = [];

allReward_R = allReward_R(relevInds);
allReward_L = allReward_L(relevInds);
allRewards = allReward_R; allRewards(~isnan(allReward_L)) = allReward_L(~isnan(allReward_L)); allRewards(allReward_L == 1) = -1;
allChoice(isnan(allChoice)) = [];
relevSessionData = sessionPupilData(relevInds);

firstFrame = 118;
lastFrame = 118;
for i = 2:length(relevInds)
    tempPupil_B = sessionPupilData(relevInds(i)-1).pupilDiameterSmooth(end-firstFrame:end);
%      tempPupil_B(sessionPupilData(relevInds(i)-1).pupilNaN(end-firstFrame:end)) = NaN;
    tempPupil_A = sessionPupilData(relevInds(i)).pupilDiameterSmooth(1:lastFrame);
%      tempPupil_A(sessionPupilData(relevInds(i)).pupilNaN(1:lastFrame)) = NaN;
    relevPupilData(i,:) = [tempPupil_B tempPupil_A];
%     relevPupilData(i,:) = [sessionPupilData(relevInds(i)-1).pupilDiameterSmooth(end-firstFrame:end) sessionPupilData(relevInds(i)).pupilDiameterSmooth(1:lastFrame)];
end
% relevPupilData(1,:) = NaN(1,1+firstFrame+lastFrame);

changeChoice = [false abs(diff(allChoice)) > 0];
time = [-firstFrame:lastFrame]/20;

%%
figure
errorfill(time, nanmean(relevPupilData(changeChoice,:)), nanstd(relevPupilData(changeChoice,:))/sqrt(sum(changeChoice)),'r');
errorfill(time, nanmean(relevPupilData(~changeChoice,:)),nanstd(relevPupilData(~changeChoice,:))/sqrt(sum(~changeChoice)),'k');
xlabel('Time (s)')
ylabel('Pupil Diameter (mm)')
legend('','Explore Trials','','Exploit Trials')

%%
figure
errorfill(time, nanmean(relevPupilData(changeChoice,:)), nanstd(relevPupilData(changeChoice,:))/sqrt(sum(changeChoice)),'r');
errorfill(time, nanmean(relevPupilData([changeChoice(2:end) false],:)),nanstd(relevPupilData([changeChoice(2:end) false],:))/sqrt(sum([changeChoice(2:end) false])),'k');
xlabel('Time (s)')
ylabel('Pupil Diameter (mm)')
legend('','Explore Trials','','1 trial before Explore Trials')
%%
figure
errorfill(time, nanmean(relevPupilData(changeChoice,:)), nanstd(relevPupilData(changeChoice,:))/sqrt(sum(changeChoice)),'r');
errorfill(time, nanmean(relevPupilData([false changeChoice(1:end-1)],:)),nanstd(relevPupilData([false changeChoice(1:end-1)],:))/sqrt(sum([false changeChoice(1:end-1)])),'k');
xlabel('Time (s)')
ylabel('Pupil Diameter (mm)')
legend('','Explore Trials','','1 trial after Explore Trials')
%%
figure
errorfill(time, nanmean(relevPupilData(allChoice == 1,:)), nanstd(relevPupilData(allChoice == 1,:))/sqrt(sum(allChoice == 1)),'r');
errorfill(time, nanmean(relevPupilData(allChoice == -1,:)), nanstd(relevPupilData(allChoice == -1,:))/sqrt(sum(allChoice == -1)),'k');

%%
figure
errorfill(time, nanmean(relevPupilData([changeChoice(2:end) false],:)),nanstd(relevPupilData([changeChoice(2:end) false],:))/sqrt(sum([changeChoice(2:end) false])),'r');
errorfill(time, nanmean(relevPupilData([false changeChoice(1:end-1)],:)),nanstd(relevPupilData([false changeChoice(1:end-1)],:))/sqrt(sum([false changeChoice(1:end-1)])),'k');
xlabel('Time (s)')
ylabel('Pupil Diameter (mm)')
legend('','1 trial before Explore Trials','','1 trial before Exploit Trials')

%%
changeChoice = [false abs(diff(allChoice)) > 0];
changeHistogram = [];
for i = find(changeChoice == 1)
    if allChoice(i) == 1 % if a right lick
        temp = 0;
        goBack = 1;
        while (i - goBack > 0) && allChoice(i-goBack) == -1 && allReward_L(i-goBack) == 0 % if previous trial was a L lick AND resulted in no reward
            temp = temp + 1;
            goBack = goBack + 1;
        end
        changeHistogram = [changeHistogram temp];
    elseif allChoice(i) == -1 %if a left lick
        temp = 0;
        goBack = 1;
        while (i - goBack > 0) && allChoice(i-goBack) == 1 && allReward_R(i-goBack) == 0 % if previous trial was a R lick AND resulted in no reward
            temp = temp + 1;
            goBack = goBack + 1;
        end
        changeHistogram = [changeHistogram temp];
    end
end

figure; histogram(changeHistogram,0:1:30)%,'Normalization','probability')

%%
changeChoiceInds = find(changeChoice);
changeChoiceInds(changeHistogram <=3 ) = [];
newChangeChoice = zeros(1,length(changeChoice));
if ~isempty(changeChoiceInds)
    newChangeChoice(changeChoiceInds) = 1;
    newChangeChoice = logical(newChangeChoice);
else
    newChangeChoice = changeChoice;
    sprintf('Regular changeChoice');
end
figure
errorfill(time, nanmean(relevPupilData(newChangeChoice,:)), nanstd(relevPupilData(newChangeChoice,:))/sqrt(sum(newChangeChoice)),'r');
errorfill(time, nanmean(relevPupilData(~newChangeChoice,:)),nanstd(relevPupilData(~newChangeChoice,:))/sqrt(sum(~newChangeChoice)),'k');
xlabel('Time (s)')
ylabel('Pupil Diameter (mm)')
legend('','Explore Trials','','Exploit Trials')

%%
figure; hold on
plot([sessionPupilData.pupilTime], [sessionPupilData.pupilDiameterSmooth],'k','linewidth',2);
for i = 1:length(sessionPupilData)
    if ~isnan(sessionPupilData(i).rewardR)
        if sessionPupilData(i).rewardR == 0 % blue
            plot([sessionPupilData(i).rewardTime sessionPupilData(i).rewardTime],[0.5 1.5], 'Color', [173,216,230]/255,'linewidth',2)
        elseif sessionPupilData(i).rewardR == 1
            plot([sessionPupilData(i).rewardTime sessionPupilData(i).rewardTime],[0.5 1.5], 'Color', [0 0 255]/255,'linewidth',2)
        end
    elseif ~isnan(sessionPupilData(i).rewardL)
        if sessionPupilData(i).rewardL == 0
            plot([sessionPupilData(i).rewardTime sessionPupilData(i).rewardTime],[0.5 1.5], 'Color', [255,165,0]/255,'linewidth',2)
        elseif sessionPupilData(i).rewardL == 1
            plot([sessionPupilData(i).rewardTime sessionPupilData(i).rewardTime],[0.5 1.5], 'Color', [128 0 0]/255,'linewidth',2)
        end
    end
end

ylim([0.6 1.0])
xlim([1157456 1157456+200000])