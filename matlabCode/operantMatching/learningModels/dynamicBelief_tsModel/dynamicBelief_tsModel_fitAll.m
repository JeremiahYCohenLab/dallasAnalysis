function model = dynamicBelief_tsModel_fitAll(filename, testFlag, shortCircuitFlag)

if nargin < 3
    shortCircuitFlag = 0;
end

if nargin < 2
    testFlag = 0;
end


[sessionData, unCorrectedBlockSwitch, out] = loadBehavioralData(filename);
behavStruct = parseBehavioralData(sessionData, unCorrectedBlockSwitch);
[root, sep] = currComputer_operantMatching();

outcome = abs([behavStruct.allReward_R; behavStruct.allReward_L])';
choice = abs([behavStruct.allChoice_R; behavStruct.allChoice_L])';

if shortCircuitFlag == 1 % break out of this function
    model.shortCircuit.bestParams = NaN;
    model.shortCircuit.LH = NaN;
    model.shortCircuit.BIC = NaN;
    model.shortCircuit.probChoice = NaN(length(choice), 2);
    return
end

% Initialize models
modelNames = {'base', 'basePrior'};
startValueCSVs = {'dynamicBelief_tsModel_base_startValues.csv', 'dynamicBelief_tsModel_basePrior_startValues.csv'};
              
% Set up optimization problem
options = optimset('Algorithm', 'interior-point','ObjectiveLimit',...
    -1.000000000e+300,'TolFun',1e-15, 'Display','off');

alpha_range = [0 1];
a_range = [1 10];
b_range = [1 10];

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
    if strcmp(modelNames{currMod}, 'base')
        b=[ alpha_range(2); 
           -alpha_range(1)];
        parfor r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@dynamicBelief_tsModel_base, startValues(r, :), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~, bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).posterior_right, ...
            model.(modelNames{currMod}).posterior_left] = ...
            dynamicBelief_tsModel_base(model.(modelNames{currMod}).bestParams, choice, outcome);
    elseif  strcmp(modelNames{currMod}, 'basePrior')
        b=[ alpha_range(2); a_range(2); b_range(2);
           -alpha_range(1); -a_range(1); -b_range(1)];
        parfor r = 1:runs
            disp(r)
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@dynamicBelief_tsModel_basePrior, startValues(r, :), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~, bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).posterior_right, ...
            model.(modelNames{currMod}).posterior_left] = ...
            dynamicBelief_tsModel_basePrior(model.(modelNames{currMod}).bestParams, choice, outcome);
    else 
        error('Dynamic belief model: Model name not found')
    end
    
    model.(modelNames{currMod}).LH = -1 * LH(bestFit, :);
    model.(modelNames{currMod}).BIC = log(length(outcome))*numParam - 2*model.(modelNames{currMod}).LH;
    
    bestHess = squeeze(hess(bestFit, :, :));
    model.(modelNames{currMod}).CIvals = sqrt(diag(inv(bestHess)))'*1.96;
    model.(modelNames{currMod}).exitFl = exitFl(bestFit, :);
end