function [mdl] = simulationLogReg(varargin)


%task and model parameters
p = inputParser;
% default parameters if none given
p.addParameter('bestParams', []);
p.addParameter('runs', 5);
p.addParameter('maxTrials', 200);
p.addParameter('randomSeed', 27);
p.addParameter('oppoFlag', 0);
p.addParameter('figFlag', 0);
p.addParameter('compareFlag', 1);
p.addParameter('animalInfo', {'goodBehDays.xlsx','CG14','preS',0});
p.parse(varargin{:});

tMax = 12;
rwdMatx = [];
noRwdMatx = [];
combinedAllChoice_R = [];
rSeed = p.Results.randomSeed;

if p.Results.compareFlag
    [actual, tMax] = combineLogReg_opMD(p.Results.animalInfo{1}, p.Results.animalInfo{2}, p.Results.animalInfo{3}, p.Results.animalInfo{4});
end

for i = 1:p.Results.runs
    rSeed = rSeed + 1;
    if p.Results.oppoFlag
        [rBar, allRewards, allChoices] = qLearningModel_opponency_sim('maxTrials', p.Results.maxTrials, 'bestParams', p.Results.bestParams,...
            'randomSeed', rSeed);
        close;
    else
        [allRewards, allChoices] = qLearningModel_sim('maxTrials', p.Results.maxTrials, 'bestParams', p.Results.bestParams,...
            'randomSeed', rSeed);
        close;
    end
    
    allNoRewards = allChoices;
    allNoRewards(find(allRewards == 1)) = 0;
    allNoRewards(find(allRewards == -1)) = 0;
    
    allChoice_R = allChoices;
    allChoice_R(find(allChoice_R == -1)) = 0;
    
    rwdMatxTmp = [];
    noRwdMatxTmp = [];
    for j = 1:tMax
        rwdMatxTmp(j,:) = [NaN(1,j) allRewards(1:end-j)];
        noRwdMatxTmp(j,:) = [NaN(1,j) allNoRewards(1:end-j)];
    end
    
   
    rwdMatx = [rwdMatx NaN(tMax,100) rwdMatxTmp];
    noRwdMatx = [noRwdMatx NaN(tMax,100) noRwdMatxTmp];
    combinedAllChoice_R = [combinedAllChoice_R NaN(1,100) allChoice_R];
end


mdl = fitglm([rwdMatx' noRwdMatx'], combinedAllChoice_R,'distribution','binomial','link','logit');
if p.Results.figFlag
    figure;
    subplot(1,2,1); hold on;
    relevInds = 2:tMax+1;
    coefVals = mdl.Coefficients.Estimate(relevInds);
    CIbands = coefCI(mdl);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color',[0.7 0 1],'linewidth',2)
    if p.Results.compareFlag
        coefVals = actual.Coefficients.Estimate(relevInds);
        CIbands = coefCI(actual);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color',[0.7 0.5 1],'linewidth',2)
        legend('simulated', 'actual')
    end
    xlabel('Reward n Trials Back')
    ylabel('\beta Coefficient')
    xlim([0.5 tMax+0.5])

    subplot(1,2,2); hold on;
    relevInds = tMax+2:length(mdl.Coefficients.Estimate);
    coefVals = mdl.Coefficients.Estimate(relevInds);
    CIbands = coefCI(mdl);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'b','linewidth',2)
    if p.Results.compareFlag
        coefVals = actual.Coefficients.Estimate(relevInds);
        CIbands = coefCI(actual);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar((1:tMax)+0.2,coefVals,errorL,errorU,'Color',[0.5 0.5 1],'linewidth',2)
        legend('simulated', 'actual')
    end
    xlabel('No Reward n Trials Back')
    ylabel('\beta Coefficient')
    xlim([0.5 tMax+0.5])
%    title(sprintf(['v = %i'], p.Results.bestParams(5)));
end
