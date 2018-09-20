function [optimalParams, LH, exitFl, hess, P, V, actions] = actorCritic(fileOrFolder)


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

alphaValue_PPE_range = [0 1];
alphaAction_range = [0 1];
alphaForgetAction_range = [0 1];
beta_range = [0 double(intmax)];
bias_range = [-5 5];
kappa_range = [0 0]; %[-5 5];
alphaValue_NPE_range = [0 1];
alphaForgetValue_range = [1 1];


startValues = [0.5, 0.5, 0.5, 10, 0, 0, 0.5, 0.5];
A=[eye(length(startValues)); -eye(length(startValues))];
b=[ alphaValue_PPE_range(2);  alphaAction_range(2);  alphaForgetAction_range(2);  beta_range(2);  bias_range(2);  kappa_range(2);  alphaValue_NPE_range(2);  alphaForgetValue_range(2);
   -alphaValue_PPE_range(1); -alphaAction_range(1); -alphaForgetAction_range(1); -beta_range(1); -bias_range(1); -kappa_range(1); -alphaValue_NPE_range(1); -alphaForgetValue_range(1);];


[optimalParams(1,:), LH(1,:), exitFl(1,:), ~, ~, ~, hess] = fmincon(@actorCriticLearning,startValues,...
            A,b,[],[],[],[],[],options,choice,outcome,sessionChangeInd);
        
hess = sqrt(diag(inv(hess)))'*1.96;
        
[P, V, actions] = plotFigure(fileOrFolder, optimalParams,choice,outcome,sessionChangeInd,behavStruct,blockSwitch);


function [LH, P, V, actions] = actorCriticLearning(startValues, choice, outcome, sessionChangeInd)

alphaValue_PPE = startValues(1);
alphaAction = startValues(2);
alphaForgetAction = startValues(3);
beta = startValues(4);
bias = startValues(5);
kappa = startValues(6);
alphaValue_NPE = startValues(7);
alphaForgetValue = startValues(8);

trials = length(choice);
V = zeros(trials,1);
actions = zeros(trials,2);

Rmin1 = zeros(trials,1);
Lmin1 = zeros(trials,1);

% Call learning rule
for t = 1 : (trials-1)
    if any(t == sessionChangeInd) % if the session has changed
        V(t) = 0; %  restart values
        actions(t,:) = 0; % restart action parameters
    end
    if choice(t, 1) == 1 % right choice
        if outcome(t,1) == 1
            V(t+1) = alphaForgetValue*V(t) + alphaValue_PPE*(outcome(t,1) - V(t)); % update global value
        else
            V(t+1) = alphaForgetValue*V(t) + alphaValue_NPE*(0 - V(t));
        end
        
        actions(t+1, 1) = alphaForgetAction*actions(t, 1) + alphaAction*(outcome(t, 1) - V(t+1));
        actions(t+1, 2) = alphaForgetAction*actions(t, 2);
        
        Rmin1(t+1) = 1;
        Lmin1(t+1) = 0;
    else % left choice
        if outcome(t,2) == 1
            V(t+1) = alphaForgetValue*V(t) + alphaValue_PPE*(outcome(t,2) - V(t));
        else
            V(t+1) = alphaForgetValue*V(t) + alphaValue_NPE*(0 - V(t));
        end

        
        actions(t+1, 2) = alphaForgetAction*actions(t, 2) + alphaAction*(outcome(t, 2) - V(t+1));
        actions(t+1, 1) = alphaForgetAction*actions(t, 1);
        
        Rmin1(t+1) = 0;
        Lmin1(t+1) = 1;
    end
end

% Call softmax rule
P = logistic([beta*(actions(:,1)-actions(:,2))+kappa*(Rmin1-Lmin1)+bias, beta*(actions(:,2)-actions(:,1))+kappa*(Lmin1-Rmin1)-bias]);

% To calculate likelihood:
Pf = 10 * choice.* P;
nonzeros = find(Pf ~= 0);
Pf = Pf(nonzeros);
% LH = -1 * prod(Pf,1);
LH = -1 * sum(log(Pf));
end


function [probChoice, V, actions] = plotFigure(fileOrFolder, optimalParams, choice, outcome, sessionChangeInd, behavStruct, blockSwitch)
    
    [~,probChoice, V, actions] = actorCriticLearning(optimalParams, choice, outcome, sessionChangeInd);

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
            'alphaValue_PPE',num2str(optimalParams(1),2), ...
            'alphaAction',num2str(optimalParams(2),2), ...
            'alphaForgetAction',num2str(optimalParams(3),2), ...
            'beta',num2str(optimalParams(4),2), ...
            'bias',num2str(optimalParams(5),2), ...
            'kappa',num2str(optimalParams(6),2)))

    for i = 1:length(blockSwitch)
        bs_loc = blockSwitch(i);
        plot([bs_loc bs_loc],[-1 1],'--','Color',[30 144 255]./255)
    end

    title(fileOrFolder)

    subplot(3,1,2); hold on
    plot(V,'linewidth',2); 
    ylabel('Value')
    xlabel('Trials')
    xlim([1 length(behavStruct.allChoice_R)])
    ylim([0 max(max(V))])

    for i = 1:length(blockSwitch)
        bs_loc = blockSwitch(i);
        plot([bs_loc bs_loc],[0 1],'--','Color',[30 144 255]./255)
    end
    subplot(3,1,3); hold on
    plot(actions,'linewidth',2)
    ylim([min(min(actions)) max(max(actions))])
end

end