function [model] = directActor_fitAll(choice, outcome, testFlag)

if nargin < 3
    testFlag = 0;
end

% Initialize models
modelNames = {'oneParam','twoParams_alphaForget','twoParams_bias','threeParams'};
startValueCSVs = {'directActorModel_1param_startValues.csv', 'directActorModel_2params_alphaForget_startValues.csv', ...
                  'directActorModel_2params_bias_startValues.csv', 'directActorModel_3params_startValues.csv'};
              
% Set up optimization problem
options = optimset('Algorithm', 'interior-point','ObjectiveLimit',...
    -1.000000000e+300,'TolFun',1e-15, 'Display','off');

alphaLearn_range = [0 double(intmax)];
alphaForget_range = [0 1];
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
    if strcmp(modelNames{currMod}, 'oneParam')
        b=[ alphaLearn_range(2);
           -alphaLearn_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@directActorModel_1param, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).actionPreferences] = directActorModel_1param(model.(modelNames{currMod}).bestParams, choice, outcome);
        
    elseif strcmp(modelNames{currMod}, 'twoParams_alphaForget')
        b=[ alphaLearn_range(2);  alphaForget_range(2);
           -alphaLearn_range(1); -alphaForget_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@directActorModel_2params_alphaForget, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).actionPreferences] = directActorModel_2params_alphaForget(model.(modelNames{currMod}).bestParams, choice, outcome);
    elseif strcmp(modelNames{currMod}, 'twoParams_bias')
        b=[ alphaLearn_range(2);  bias_range(2);
           -alphaLearn_range(1); -bias_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@directActorModel_2params_bias, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).actionPreferences] = directActorModel_2params_bias(model.(modelNames{currMod}).bestParams, choice, outcome);
    elseif strcmp(modelNames{currMod}, 'threeParams')
        b=[ alphaLearn_range(2);  alphaForget_range(2);  bias_range(2);
           -alphaLearn_range(1); -alphaForget_range(1); -bias_range(1)];
        for r = 1:runs
            [allParams(r, :), LH(r, :), exitFl(r, :), ~, ~, ~, hess(r, :, :)] = ...
                fmincon(@directActorModel_3params, startValues(r,:), A, b, [], [], [], [], [], options, choice, outcome);
        end
        [~,bestFit] = min(LH);
        model.(modelNames{currMod}).bestParams = allParams(bestFit, :);
        [~, model.(modelNames{currMod}).probChoice, model.(modelNames{currMod}).actionPreferences] = directActorModel_3params(model.(modelNames{currMod}).bestParams, choice, outcome);
    else
        error('Direct actor: Model name not found')
    end
 
    model.(modelNames{currMod}).LH = -1 * LH(bestFit, :);
    model.(modelNames{currMod}).BIC = log(length(outcome))*numParam - 2*model.(modelNames{currMod}).LH;
    
    bestHess = squeeze(hess(bestFit, :, :));
    model.(modelNames{currMod}).CIvals = sqrt(diag(inv(bestHess)))'*1.96;
    model.(modelNames{currMod}).exitFl = exitFl(bestFit, :);
end