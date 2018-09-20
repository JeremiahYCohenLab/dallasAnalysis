function [optimalParams, LH, exitFl, hess] = directActor_withValue(fileOrFolder)


% Pull out reward and choice data
if ischar(fileOrFolder)
    fileOrFolder = {fileOrFolder};
end
outcome = [];
choice = [];
sessionChangeInd = [];
for i = 1:length(fileOrFolder)    
    [sessionData, blockSwitch] = loadBehavioralData(fileOrFolder{i});
    behavStruct = parseBehavioralData(sessionData, blockSwitch);

    outcome = [outcome; abs([behavStruct.allReward_R; behavStruct.allReward_L])'];
    choice = [choice; abs([behavStruct.allChoice_R; behavStruct.allChoice_L])'];
    sessionChangeInd = [sessionChangeInd length(outcome)];
end

% Set up optimization problem
options = optimset('Algorithm', 'interior-point','ObjectiveLimit',...
    -1.000000000e+300,'TolFun',1e-15, 'Display','off');

alphaAction_range = [0 1];
alphaForgetAction_range = [0 1];
beta_range = [0 double(intmax)];
bias_range = [-5 5];
kappa_range = [0 0]; %[-5 5];
alphaValue_range = [0.2 0.2];


startValues = [0.5, 0.5, 10, 0, 0, 0.1];
A=[eye(length(startValues)); -eye(length(startValues))];
b=[ alphaAction_range(2);  alphaForgetAction_range(2);  beta_range(2);  bias_range(2);  kappa_range(2);  alphaValue_range(2);
   -alphaAction_range(1); -alphaForgetAction_range(1); -beta_range(1); -bias_range(1); -kappa_range(1); -alphaValue_range(1)];


[optimalParams(1,:), LH(1,:), exitFl(1,:), ~, ~, ~, hess] = fmincon(@directActorLearning,startValues,...
            A,b,[],[],[],[],[],options,choice,outcome,sessionChangeInd);
        

plotFigure(fileOrFolder, optimalParams,choice,outcome,sessionChangeInd,behavStruct,blockSwitch);


function [LH, P, V, actions] = directActorLearning(startValues, choice, outcome, sessionChangeInd)

    
alphaAction = startValues(1);
alphaForgetAction = startValues(2);
beta = startValues(3);
bias = startValues(4);
kappa = startValues(5);
alphaValue = startValues(6);

trials = length(choice);
actions = zeros(trials,2);
P = zeros(trials,2);
V = zeros(trials,1);

Rmin1 = zeros(trials,1);
Lmin1 = zeros(trials,1);

% Call learning rule
t = 1;
P(t,:) = logistic([beta*(actions(t,1) - actions(t,2)) + kappa*(Rmin1(t)-Lmin1(t)) + bias, ...
                   beta*(actions(t,2) - actions(t,1)) + kappa*(Lmin1(t)-Rmin1(t)) - bias]);
for t = 1 : (trials-1)
    if any(t == sessionChangeInd) % if the session has changed
        V(t) = 0;
        actions(t,:) = 0; % restart action parameters
    end
    
    if choice(t, 1) == 1 % right choice
        V(t+1) = V(t) + alphaValue*(outcome(t,1) - V(t)); % update global value
        
        actions(t+1, 1) = alphaForgetAction*actions(t, 1) + alphaAction*(1 - P(t,1))*(outcome(t,1) - V(t+1));
        actions(t+1, 2) = alphaForgetAction*actions(t, 2) - alphaAction*P(t,2)*(outcome(t,1) - V(t+1));

        
        Rmin1(t+1) = 1;
        Lmin1(t+1) = 0;
    else % left choice
       V(t+1) = V(t) + alphaValue*(outcome(t,2) - V(t)); % update global value
       
        actions(t+1, 2) = alphaForgetAction*actions(t, 2) + alphaAction*(1 - P(t,2))*(outcome(t,2) - V(t+1));
        actions(t+1, 1) = alphaForgetAction*actions(t, 1) - alphaAction*P(t,1)*(outcome(t,2) - V(t+1));

        
        Rmin1(t+1) = 0;
        Lmin1(t+1) = 1;
    end
    P(t+1,:) = logistic([beta*(actions(t+1,1) - actions(t+1,2)) + kappa*(Rmin1(t) - Lmin1(t)) + bias, ...
                         beta*(actions(t+1,2) - actions(t+1,1)) + kappa*(Lmin1(t) - Rmin1(t)) - bias]);
end

% To calculate likelihood:
Pf = 10 * choice.* P;
nonzeros = find(Pf ~= 0);
Pf = Pf(nonzeros);
% LH = -1 * prod(Pf,1);
LH = -1 * sum(log(Pf));
end


function plotFigure(fileOrFolder, optimalParams, choice, outcome, sessionChangeInd, behavStruct, blockSwitch)
    
    [~,probChoice, V, actions] = directActorLearning(optimalParams, choice, outcome, sessionChangeInd);

    normKern = normpdf(-15:15,0,2);
    normKern = normKern / sum(normKern);
    xVals = (1:(length(normKern) + length(behavStruct.allChoices) - 1)) - round(length(normKern)/2);

    figure; 
    subplot(3,1,1); hold on
    plot(xVals, conv(behavStruct.allChoices,normKern),'k','linewidth',2);
    plot(xVals, conv(behavStruct.allRewards,normKern),'--','Color',[100 100 100]./255,'linewidth',2)
    xlabel('Trials')
    ylabel('<-- Left       Right -->')
    xlim([1 length(behavStruct.allChoice_R)])
    ylim([-1 1])

    plot(probChoice(:,1)*2-1,'-','Color',[148,0,211]./255,'linewidth',1)

    legend('Choices','Rewards', ...
            sprintf('%s\n%s: %s\n%s: %s\n%s: %s\n%s: %s\n%s: %s\n%s: %s', 'Q-learning Model',....
            'alphaAction',num2str(optimalParams(1),2), ...
            'alphaForgetAction',num2str(optimalParams(2),2), ...
            'beta',num2str(optimalParams(3),2), ...
            'bias',num2str(optimalParams(4),2), ...
            'kappa',num2str(optimalParams(5),2), ...
            'alphaValue',num2str(optimalParams(6),2)))

    for i = 1:length(blockSwitch)
        bs_loc = blockSwitch(i);
        plot([bs_loc bs_loc],[-1 1],'--','Color',[30 144 255]./255)
    end

    title(fileOrFolder)

    subplot(3,1,2); hold on
    plot(V,'linewidth',2); 
    ylabel('Values')
    xlabel('Trials')
    xlim([1 length(behavStruct.allChoice_R)])
    ylim([min(min(V)) max(max(V))])

    for i = 1:length(blockSwitch)
        bs_loc = blockSwitch(i);
        plot([bs_loc bs_loc],[0 1],'--','Color',[30 144 255]./255)
    end
    subplot(3,1,3); hold on
    plot(actions,'linewidth',2);
    legend('Right','Left')
end

end