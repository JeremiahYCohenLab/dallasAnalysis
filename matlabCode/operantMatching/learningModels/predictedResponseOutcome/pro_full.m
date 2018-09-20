function [LH, probChoice, V, actions] = pro_full(startValues, choice, outcome)

alphaValue = startValues(1);
alphaLearn_action = startValues(2);
alphaForget_action = startValues(3);
beta = startValues(4);

trials = length(choice);
V = zeros(trials,1);
actions = zeros(trials,2);

% Call learning rule
for t = 1 : (trials-1)
    if choice(t, 1) == 1 % right choice
        errorV = outcome(t, 1) - V(t);
        alphaValue_corrected = alphaValue / (1 + abs(errorV));        
        V(t+1) = V(t) + alphaValue_corrected*errorV; % update global value
        
        errorA = outcome(t, 1) - V(t+1);
        alphaLearn_action_corrected = alphaLearn_action / (1 + abs(errorA));
        actions(t+1, 1) = alphaForget_action*actions(t, 1) + alphaLearn_action_corrected*errorA;
        actions(t+1, 2) = alphaForget_action*actions(t, 2);
        
    else % left choice
        errorV = outcome(t, 2) - V(t);
        alphaValue_corrected = alphaValue / (1 + abs(errorV));
        V(t+1) = V(t) + alphaValue_corrected / (1 + abs(errorV));
        
        errorA = outcome(t, 2) - V(t+1);
        alphaLearn_action_corrected = alphaLearn_action / (1 + abs(errorA));
        actions(t+1, 2) = alphaForget_action*actions(t, 2) + alphaLearn_action_corrected*errorA;
        actions(t+1, 1) = alphaForget_action*actions(t, 1);

    end
end

% Call softmax rule
probChoice = logistic([beta*(actions(:,1) - actions(:,2)), ...
                       beta*(actions(:,2) - actions(:,1))]);

% To calculate likelihood:
LH = likelihood(choice, probChoice);
end