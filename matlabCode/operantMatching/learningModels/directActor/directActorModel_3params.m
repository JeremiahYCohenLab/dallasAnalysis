function [LH, probChoice, actionPreference] = directActorModel_3params(startValues, choice, outcome)

alphaLearn = startValues(1);
alphaForget = startValues(2);
bias = startValues(3);

trials = length(choice);
actionPreference = zeros(trials,2);
probChoice = zeros(trials,2);

% Call learning rule
% initialize the first probChoice
probChoice(1,:) = logistic([(actionPreference(1, 1) - actionPreference(1, 2)) + bias, ...
                            (actionPreference(1, 2) - actionPreference(1, 1)) - bias]);
for t = 1 : (trials-1)
    
    if choice(t, 1) == 1 % right choice
        actionPreference(t+1, 1) = alphaForget*actionPreference(t, 1) + alphaLearn*(1 - probChoice(t, 1))*outcome(t, 1);
        actionPreference(t+1, 2) = alphaForget*actionPreference(t, 2) - alphaLearn*probChoice(t, 2)*outcome(t, 1);
    else % left choice 
        actionPreference(t+1, 2) = alphaForget*actionPreference(t, 2) + alphaLearn*(1 - probChoice(t, 2))*outcome(t, 2);
        actionPreference(t+1, 1) = alphaForget*actionPreference(t, 1) - alphaLearn*probChoice(t, 1)*outcome(t, 2);
    end
    probChoice(t+1, :) = logistic([(actionPreference(t+1, 1) - actionPreference(t+1,2)) + bias, ...
                                   (actionPreference(t+1, 2) - actionPreference(t+1,1)) - bias]);
end

% To calculate likelihood:
LH = likelihood(choice,probChoice);
end