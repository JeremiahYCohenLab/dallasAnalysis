function [paramStruct] = qLearning_fitOpponencyFiveAll(xlFile, animal, category)

testFlag = 0; 

[root, sep] = currComputer();

[weights, dayList, ~] = xlsread(xlFile, animal);
[~,col] = find(~cellfun(@isempty,strfind(dayList, category)) == 1);
dayList = dayList(2:end,col);
endInd = find(cellfun(@isempty,dayList),1);
if ~isempty(endInd)
    dayList = dayList(1:endInd-1,:);
end

allBIC = [];
allAlphaPPE = [];
allAlphaNPE = [];
allAlphaForget = [];
allBeta = [];
allV = [];
allRbar = [];

for i = 1: length(dayList)
    filename = [dayList{i} '.asc'];

    [behSessionData, unCorrectedBlockSwitch, out] = loadBehavioralData(filename);
    behavStruct = parseBehavioralData(behSessionData, unCorrectedBlockSwitch);
    [root, sep] = currComputer_operantMatching();

    outcome = abs([behavStruct.allReward_R; behavStruct.allReward_L])';
    choice = abs([behavStruct.allChoice_R; behavStruct.allChoice_L])';


    % Initialize models
    modelNames = {'fiveParams_opponency'};
    startValueCSVs = {'qLearningModel_5params_opponency_startValues.csv'};

    % Set up optimization problem
    options = optimset('Algorithm', 'interior-point','ObjectiveLimit',...
        -1.000000000e+300,'TolFun',1e-15, 'Display','off');

    alpha_range = [0 1];
    alphaNPE_range = [0 1];
    alphaPPE_range = [0 1];
    alphaForget_range = [0 1];
    beta_range = [0 double(intmax)];
    v_range = [0 1];

    for currMod = 1:length(modelNames)
        startValues = csvread(startValueCSVs{currMod});
        if testFlag == 1
            startValues = startValues(1, :);
        end

        % initialize output variables
        runs = size(startValues, 1);
        allParams = zeros(size(startValues, 1), size(startValues, 2));
        LH = zeros(size(startValues, 1), 1);
        exitFl = zeros(size(startValues, 1), 1);
        hess = zeros(size(startValues, 1), size(startValues, 2), size(startValues, 2));
        numParam = size(startValues, 2);
        A=[eye(size(startValues, 2)); -eye(size(startValues, 2))];


        if strcmp(modelNames{currMod}, 'fiveParams_opponency')
            b=[ alphaNPE_range(2);  alphaPPE_range(2);  alphaForget_range(2); beta_range(2); v_range(2);
               -alphaNPE_range(1); -alphaPPE_range(1); -alphaForget_range(1); -beta_range(1); -v_range(1)];
            for r = 1:runs
                [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                    fmincon(@qLearningModel_5params_opponency, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
            end
            [~,bestFit] = min(LH);
            model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
            [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).Q, model.(modelNames{currMod}).pe, model.(modelNames{currMod}).rBar] = ...
                qLearningModel_5params_opponency(model.(modelNames{currMod}).bestParams, choice, outcome);
        end

        model.(modelNames{currMod}).LH = -1 * LH(bestFit, :);
        model.(modelNames{currMod}).BIC = log(length(outcome))*numParam - 2*model.(modelNames{currMod}).LH;

        bestHess = squeeze(hess(bestFit, :, :));
        model.(modelNames{currMod}).CIvals = sqrt(diag(inv(bestHess)))'*1.96;
        model.(modelNames{currMod}).exitFl = exitFl(bestFit, :);
    end
   
    
    allAlphaPPE = [allAlphaPPE model.fiveParams_opponency.bestParams(1)];
    allAlphaNPE = [allAlphaNPE model.fiveParams_opponency.bestParams(2)];
    allAlphaForget = [allAlphaForget model.fiveParams_opponency.bestParams(3)];
    allBeta = [allBeta model.fiveParams_opponency.bestParams(4)];
    allV = [allV model.fiveParams_opponency.bestParams(5)];
    
    allBIC = [allBIC model.fiveParams_opponency.BIC];
    allRbar = [allRbar model.fiveParams_opponency.rBar'];
end

paramStruct.allBIC = allBIC;
paramStruct.allAlphaPPE = allAlphaPPE;
paramStruct.allAlphaNPE = allAlphaNPE;
paramStruct.allAlphaForget = allAlphaForget;
paramStruct.allBeta = allBeta;
paramStruct.allV = allV;
paramStruct.allRbar = allRbar;
    
    