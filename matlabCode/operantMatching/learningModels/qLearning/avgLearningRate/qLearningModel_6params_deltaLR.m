function [LH, probChoice, Q, pe, rBar, alphaNPE, alphaPPE] = qLearningModel_6params_deltaLR(startValues, choice, outcome)

alphaNPEstart = startValues(1);
alphaPPEstart = startValues(2);
alphaForget = startValues(3);
beta = startValues(4);
v = startValues(5);
w = startValues(6);

trials = length(choice);
Q = zeros(trials,2);
rBar = [0; zeros(trials-1,1)];
alphaNPE = [alphaNPEstart; zeros(trials-1,1)];
alphaPPE = [alphaPPEstart; zeros(trials-1,1)];

% Call learning rule
for t = 1 : (trials-1)
    if choice(t, 1) == 1 % right choice
        Q(t+1, 2) = alphaForget*Q(t, 2);
        pe(t) = outcome(t, 1) - Q(t, 1);
        if pe(t) < 0
            Q(t+1, 1) = alphaForget*Q(t, 1) + alphaNPE(t) * pe(t);
        else
            Q(t+1, 1) = alphaForget*Q(t, 1) + alphaPPE(t) * pe(t);
        end
        rBar(t+1) = v*outcome(t, 1) + (1-v)*rBar(t);
    else % left choice
        Q(t+1, 1) = alphaForget*Q(t, 1);
        pe(t) = outcome(t, 2) - Q(t, 2);
        if pe(t) < 0
            Q(t+1, 2) = alphaForget*Q(t, 2) + alphaNPE(t) * pe(t);
        else
            Q(t+1, 2) = alphaForget*Q(t, 2) + alphaPPE(t) * pe(t);
        end
        rBar(t+1) = v*outcome(t, 2) + (1-v)*rBar(t);
    end
    alphaNPE(t+1) = alphaNPEstart + w*rBar(t+1);
    alphaPPE(t+1) = alphaPPEstart - w*rBar(t+1);
    if alphaPPE(t+1) < 0 
        alphaPPE(t+1) = 0;
    end
end

if choice(t, 1) == 1
    pe(trials) = outcome(end, 1) - Q(end, 1);
else
    pe(trials) = outcome(end, 2) - Q(end, 2);
end


% Call softmax  rule

probChoice = logistic([beta*(Q(:, 1)-Q(:, 2)), ...
                       beta*(Q(:, 2)-Q(:, 1))]);

% To calculate likelihood:
LH = likelihood(choice,probChoice);
end