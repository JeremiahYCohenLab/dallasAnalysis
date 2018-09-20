function model = stackedProbability_fitAll(choice, outcome, testFlag)

% model from Huh, Jo, Kim, Sul, Jung (2009)

if nargin < 3
    testFlag = 0;
end

% Initialize models
modelNames = {'base','baseForget','baseBias','full'};
startValueCSVs = {'spModel_base_startValues.csv','spModel_baseForget_startValues.csv', ...
                  'spModel_baseBias_startValues.csv', 'spModel_full_startValues.csv'};
              
% Set up optimization problem
options = optimset('Algorithm', 'interior-point','ObjectiveLimit',...
    -1.000000000e+300,'TolFun',1e-15, 'Display','off');

alphaLearn_range = [0 1];
alphaForget_range = [0 1];
beta_range = [0 double(intmax)];
bias_range = [-5 5];

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
        b = [ alphaLearn_range(2);  beta_range(2);
             -alphaLearn_range(1); -beta_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@spModel_base, startValues(r, :), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~, bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).X, model.(modelNames{currMod}).A, model.(modelNames{currMod}).S] = ...
            spModel_base(model.(modelNames{currMod}).bestParams, choice, outcome);
    elseif strcmp(modelNames{currMod}, 'baseForget')
        b=[ alphaLearn_range(2);  alphaForget_range(2);  beta_range(2);
           -alphaLearn_range(1); -alphaForget_range(1); -beta_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@spModel_baseForget, startValues(r, :), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~, bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).X, model.(modelNames{currMod}).A, model.(modelNames{currMod}).S] = ...
            spModel_baseForget(model.(modelNames{currMod}).bestParams, choice, outcome);
    elseif strcmp(modelNames{currMod}, 'baseBias')
        b=[ alphaLearn_range(2);  beta_range(2);  bias_range(2); 
           -alphaLearn_range(1); -beta_range(1); -bias_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@spModel_baseBias, startValues(r, :), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~, bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).X, model.(modelNames{currMod}).A, model.(modelNames{currMod}).S] = ...
            spModel_baseBias(model.(modelNames{currMod}).bestParams, choice, outcome);
    elseif strcmp(modelNames{currMod}, 'full')
        b=[ alphaLearn_range(2);  alphaForget_range(2);  beta_range(2);  bias_range(2); 
           -alphaLearn_range(1); -alphaForget_range(1); -beta_range(1); -bias_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@spModel_full, startValues(r, :), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~, bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).X, model.(modelNames{currMod}).A, model.(modelNames{currMod}).S] = ...
            spModel_full(model.(modelNames{currMod}).bestParams, choice, outcome);
    else 
        error('Stacked probability: Model name not found')
    end
    
    model.(modelNames{currMod}).LH = -1 * LH(bestFit, :);
    model.(modelNames{currMod}).BIC = log(length(outcome))*numParam - 2*model.(modelNames{currMod}).LH;
    
    bestHess = squeeze(hess(bestFit, :, :));
    model.(modelNames{currMod}).CIvals = sqrt(diag(inv(bestHess)))'*1.96;
    model.(modelNames{currMod}).exitFl = exitFl(bestFit, :);
end