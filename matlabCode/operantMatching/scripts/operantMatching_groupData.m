filename = {'mBB040d20160911','mBB040d20160912','mBB040d20160913', ...
    'mBB040d20160914','mBB040d20160915','mBB040d20160916', 'mBB040d20160917b', ...
    'mBB040d20160918','mBB040d20160919'};
allR = []; allC = [];
for i=1:length(filename)
[tmpR, tmpC] = analyzeBehavioralData_operantMatching(filename{i});
allR = [allR NaN(1,1000) tmpR];
allC = [allC NaN(1,1000) tmpC];
close all
end

%% Models
tMax = 30;
rwdMatx = [];
for i = 1:tMax
rwdMatx(i,:) = [NaN(1,i) allR(1:end-i)];
end
choiceMatx = [];
for i = 1:tMax
choiceMatx(i,:) = [NaN(1,i) allC(1:end-i)];
end
lm_rwd = fitlm([rwdMatx]', allC,'linear');
lm_choice = fitlm([choiceMatx]', allC,'linear');
lm_rwdXchoice = fitlm([rwdMatx; choiceMatx]', allC,'linear');

%% Reward
figure; subplot(2,1,1)

currModel = lm_rwdXchoice;
coefVals = currModel.Coefficients.Estimate(2:tMax+1)';
allCI = currModel.coefCI;
errorTerms = [abs(currModel.Coefficients.Estimate(2:tMax+1) - allCI(2:tMax+1,1)), abs(currModel.Coefficients.Estimate(2:tMax+1)-allCI(2:tMax+1,2))]';
errorfill(1:tMax,coefVals,errorTerms,'r');

currModel = lm_rwd;
coefVals = currModel.Coefficients.Estimate(2:tMax+1)';
allCI = currModel.coefCI;
errorTerms = [abs(currModel.Coefficients.Estimate(2:tMax+1) - allCI(2:tMax+1,1)), abs(currModel.Coefficients.Estimate(2:tMax+1)-allCI(2:tMax+1,2))]';
errorfill(1:tMax,coefVals,errorTerms,'k');

hold on; xlim([1 tMax])
plot([0 50],[0 0],'--','Color',[120 120 120]./255)

xlabel('Reward n Trials Back')
legend('Reward + Choice','','Reward Only','')
%% Choice
subplot(2,1,2)

currModel = lm_rwdXchoice;
coefVals = currModel.Coefficients.Estimate(tMax+2:end)';
allCI = currModel.coefCI;
errorTerms = [abs(currModel.Coefficients.Estimate(tMax+2:end) - allCI(tMax+2:end,1)), abs(currModel.Coefficients.Estimate(tMax+2:end)-allCI(tMax+2:end,2))]';
errorfill(1:tMax,coefVals,errorTerms,'r');

currModel = lm_choice;
coefVals = currModel.Coefficients.Estimate(2:tMax+1)';
allCI = currModel.coefCI;
errorTerms = [abs(currModel.Coefficients.Estimate(2:tMax+1) - allCI(2:tMax+1,1)), abs(currModel.Coefficients.Estimate(2:tMax+1)-allCI(2:tMax+1,2))]';
errorfill(1:tMax,coefVals,errorTerms,'k');

hold on; xlim([1 tMax])
plot([0 50],[0 0],'--','Color',[120 120 120]./255)

xlabel('Choice n Trials Back')
legend('Reward + Choice','','Choice Only','')