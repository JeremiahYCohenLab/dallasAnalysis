function [rBarTbl, rBarTbl_zero] = rBarAnalysis_opMD(xlFile, animal, category, revForFlag, plotFlag)

if nargin < 5
    plotFlag = 0;
end
if nargin < 4
    revForFlag = 0;
end

[root, sep] = currComputer();

[weights, dayList, ~] = xlsread(xlFile, animal);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

smoothKern = ones(1,10);
purp = [0.5 0 0.8];
blue = [0 0 1];


for i = 1:length(dayList)
    mdl = qLearning_fitOpponency([dayList{i} '.asc'], revForFlag);
    alphaNPE(i) = mdl.fiveParams_opponency.bestParams(1);
    alphaPPE(i) = mdl.fiveParams_opponency.bestParams(2);
    alphaForget(i) = mdl.fiveParams_opponency.bestParams(3);
    beta(i) = mdl.fiveParams_opponency.bestParams(4);
    v(i) = mdl.fiveParams_opponency.bestParams(5);
    alphaNPE_zero(i) = mdl.fiveParams_opponency_zero.bestParams(1);
    alphaPPE_zero(i) = mdl.fiveParams_opponency_zero.bestParams(2);
    alphaForget_zero(i) = mdl.fiveParams_opponency_zero.bestParams(3);
    beta_zero(i) = mdl.fiveParams_opponency_zero.bestParams(4);
    v_zero(i) = mdl.fiveParams_opponency_zero.bestParams(5);
    rwdArray = zeros(length(mdl.outcome), 1);
    rwdArray(find(mdl.outcome(:,1) == 1)) = 1; rwdArray(find(mdl.outcome(:,2) == 1)) = 1;
    rwdRate(i) = sum(rwdArray) / length(rwdArray);
    rwdSmooth = conv(rwdArray', smoothKern); rwdSmooth = rwdSmooth(1:end-(length(smoothKern)-1));
    rwdVar(i) = std(rwdSmooth);
    
    if plotFlag
        figure;
        set(gcf,'defaultAxesColorOrder',[purp; blue]);
        yyaxis left; plot(mdl.fiveParams_opponency.rBar, '-', 'Color', purp, 'LineWidth', 2);
        ylabel('rBar')
        yyaxis right; plot(rwdSmooth,  '-', 'Color', blue, 'LineWidth', 2);
        ylabel('smoothed rewards')
        xlabel('trials')
        suptitle(dayList{i})
    end
end

rBarTbl = table(rwdRate', rwdVar', alphaNPE', alphaPPE', alphaForget', beta', v',... 
    'VariableNames', {'rwdRate' 'rwdVar' 'alphaNPE' 'alphaPPE' 'alphaForget' 'beta' 'v'}, 'RowNames', [dayList]);
rBarTbl_zero = table(rwdRate', rwdVar', alphaNPE_zero', alphaPPE_zero', alphaForget_zero', beta_zero', v_zero',... 
    'VariableNames', {'rwdRate' 'rwdVar' 'alphaNPE' 'alphaPPE' 'alphaForget' 'beta' 'v'}, 'RowNames', [dayList])


if plotFlag
    figure;
    purp = [0.5 0 0.8];
    blue = [0 0 1];
    set(gcf,'defaultAxesColorOrder',[purp; blue]);

    subplot(3,2,1);
    yyaxis left; plot(rBarTbl_zero.rwdRate, '-', 'Color', purp);
    ylabel('reward rate')
    yyaxis right; plot(log(rBarTbl_zero.v), '-', 'Color', blue);
    ylabel('ln(v)')
    xlabel('session')
    title('session reward rate')

    subplot(3,2,3); 
    crosscorr(rBarTbl_zero.rwdRate,log(rBarTbl_zero.v))

    subplot(3,2,5); 
    binEdges = [min(rBarTbl_zero.rwdRate) : (max(rBarTbl_zero.rwdRate) - min(rBarTbl_zero.rwdRate)) / 7 : max(rBarTbl_zero.rwdRate)];
    for j = 1:length(binEdges)-1
        vTmp = rBarTbl_zero.v(find(rBarTbl_zero.rwdRate > binEdges(j) & rBarTbl_zero.rwdRate < binEdges(j+1)));
        meanV(j) = mean(log(vTmp));
        semV(j) = std(log(vTmp)) / sqrt(length(vTmp));
    end
    xVals = binEdges(1:end-1) + diff(binEdges)/2; 
    errorbar(xVals, meanV, semV);
    xlabel('reward rate')
    ylabel('ln(v)')

    subplot(3,2,2);
    yyaxis left; plot(rBarTbl_zero.rwdVar, '-', 'Color', purp);
    ylabel('reward rate')
    yyaxis right; plot(log(rBarTbl_zero.v), '-', 'Color', blue);
    ylabel('ln(v)')
    xlabel('session')
    title('session reward rate variability')

    subplot(3,2,4); 
    crosscorr(rBarTbl_zero.rwdVar,log(rBarTbl_zero.v))

    subplot(3,2,6); 
    binEdges = [min(rBarTbl_zero.rwdVar) : (max(rBarTbl_zero.rwdVar) - min(rBarTbl_zero.rwdVar)) / 7 : max(rBarTbl_zero.rwdVar)];
    for j = 1:length(binEdges)-1
        vTmp = rBarTbl_zero.v(find(rBarTbl_zero.rwdVar > binEdges(j) & rBarTbl_zero.rwdVar < binEdges(j+1)));
        meanV(j) = mean(log(vTmp));
        semV(j) = std(log(vTmp)) / sqrt(length(vTmp));
    end
    xVals = binEdges(1:end-1) + diff(binEdges)/2; 
    errorbar(xVals, meanV, semV);
    xlabel('reward rate variability')
    ylabel('ln(v)')
end