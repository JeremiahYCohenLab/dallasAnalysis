function [LH, probChoice, posterior_right, posterior_left] = dynamicBelief_tsModel_base(startValues, choice, outcome)

alpha = startValues(1); % uncertainty parameter

a = 1;
b = 1;
gamma = linspace(0, 1, 1e3);
p0 = gamma.^(a - 1).*(1 - gamma).^(b - 1); p0 = p0/sum(p0);
iterposterior_Left =  p0; iterposterior_Left = iterposterior_Left/sum(iterposterior_Left);
iterposterior_Right = p0; iterposterior_Right = iterposterior_Right/sum(iterposterior_Right);

trials = length(choice);

posterior_right = NaN(trials, length(gamma));
posterior_left = NaN(trials, length(gamma));
posterior_right(1, :) = iterposterior_Right;
posterior_left(1, :) = iterposterior_Left;

temp_pR = 0;
for i = 1:length(gamma)
    temp_pR = temp_pR + sum(iterposterior_Left(1:i))*iterposterior_Right(i);
end
probChoice(1, :) = [temp_pR, (1 - temp_pR)];

% Call learning rule
for t = 1 : (trials-1)
    r = sum(outcome(t, :)); % 1 or 0 for outcome
    lhood = (gamma.^r).*(1 - gamma).^(1 - r);

    if choice(t, 1) == 1 % right choice
        iterposterior_Right = lhood.*iterposterior_Right;
        iterposterior_Right = iterposterior_Right/sum(iterposterior_Right);
    elseif choice(t, 2) == 1 % left choice
        iterposterior_Left = lhood.*iterposterior_Left;
        iterposterior_Left = iterposterior_Left/sum(iterposterior_Left);
    end
    
    iterposterior_Left = alpha*iterposterior_Left + (1 - alpha)*p0;
    iterposterior_Left = iterposterior_Left/sum(iterposterior_Left);
    
    iterposterior_Right = alpha*iterposterior_Right + (1 - alpha)*p0;
    iterposterior_Right = iterposterior_Right/sum(iterposterior_Right);
    
    % Calculate P(Right > Left); Thompson sampling
    temp_pR = 0;
    for i = 1:length(gamma)
        temp_pR = temp_pR + sum(iterposterior_Left(1:i))*iterposterior_Right(i);
    end
    probChoice(t + 1, :) = [temp_pR, (1 - temp_pR)];
    
    posterior_left(t + 1, :) = iterposterior_Left;
    posterior_right(t + 1, :) = iterposterior_Right;
     
end

% To calculate likelihood:
LH = likelihood_opMD(choice, probChoice);
end

