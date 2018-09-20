function plotLickRegressionModel(xlFile, animals, categories, varargin)

p = inputParser;
% default parameters if none given
p.addParameter('revForFlag',0)
p.addParameter('makeFigFlag', 0)
p.parse(varargin{:});

timeMax = 76600;
binSize = 6300;
timeBinEdges = [1000:binSize:timeMax];
tMax = length(timeBinEdges) - 1;
relevInds = 2:tMax+1;

figure; hold on;
for i = 1:length(animals)
 subplot(1,2,i); hold on;
 mdlTmp = combineLogRegLickLatTime_opMD(xlFile, animals{i}, categories{1}, 'revForFlag',...
     p.Results.revForFlag, 'makeFigFlag', p.Results.makeFigFlag);
 
    coefVals = mdlTmp.Coefficients.Estimate(relevInds);
    CIbands = coefCI(mdlTmp);
    errorL = abs(coefVals - CIbands(relevInds,1));
    errorU = abs(coefVals - CIbands(relevInds,2));
    errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'-c','linewidth',2)

    if length(categories) > 1
        mdlTmp = combineLogRegLickLatTime_opMD(xlFile, animals{i}, categories{2}, 'revForFlag',...
            p.Results.revForFlag, 'makeFigFlag', p.Results.makeFigFlag);
        coefVals = mdlTmp.Coefficients.Estimate(relevInds);
        CIbands = coefCI(mdlTmp);
        errorL = abs(coefVals - CIbands(relevInds,1));
        errorU = abs(coefVals - CIbands(relevInds,2));
        errorbar(((1:tMax)*binSize/1000),coefVals,errorL,errorU,'Color', [0.7 0 1],'linewidth',2)
        legend('pre-lesion', 'post-lesion')
    end

    xlabel('reward n seconds back')
    ylabel('\beta Coefficient')
    xlim([0 (tMax*binSize/1000 + 5)])

    title(animals{i})
end
    
suptitle('LRM: rewards on lick latency')