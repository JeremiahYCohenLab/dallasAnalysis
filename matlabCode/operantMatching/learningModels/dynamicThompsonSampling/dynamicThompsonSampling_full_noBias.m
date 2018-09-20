function [LH, probChoice, betaDist_right, betaDist_left] = dynamicThompsonSampling_full_noBias(startValues, choice, outcome)

alphaLearn_reward = startValues(1); % scales how much to add to "a" parameter of beta dist
alphaLearn_noreward = startValues(2); % scales how much to add to "b" parameter
alphaForget = startValues(3); % decays beta dist to priors
priorA = startValues(4); % sets prior "a"
priorB = startValues(5); % sets prior "b"

trials = length(choice);

betaDist_right = ones(trials,2);
betaDist_left = ones(trials,2);

betaDist_right(1,:) = [priorA priorB];
betaDist_left(1,:) = [priorA priorB];


% % moment match beta distribution to normal distribution
% uR = betaDist_right(1,1) / (betaDist_right(1,1)+ betaDist_right(1,2));
% uL = betaDist_left(1,1)  / (betaDist_left(1,1) + betaDist_left(1,2));
% sigR = (betaDist_right(1,1) * betaDist_right(1,2)) / ((betaDist_right(1,1) + betaDist_right(1,2))^2 * ...
%         (betaDist_right(1,1) + betaDist_right(1,2) + 1));
% sigL = (betaDist_left(1,1) * betaDist_left(1,2)) / ((betaDist_left(1,1) + betaDist_left(1,2))^2 * ...
%         (betaDist_left(1,1) + betaDist_left(1,2) + 1));
% Z = ((uR - uL) / ((sigR + sigL)^(1/2)));
% probChoice(1,:) = [normcdf(Z) (1 - normcdf(Z))];

Ra = betaDist_right(1, 1);
Rb = betaDist_right(1, 2);
La = betaDist_left(1, 1);
Lb = betaDist_left(1, 2);
fun = @(x)(x.^(Ra - 1).*(1 - x).^(Rb - 1))/beta(Ra, Rb).*betacdf(x, La, Lb);
tP = integral(fun, 0, 1); % probability of right choice

probChoice(1, :) = [tP, (1-tP)];

% Call learning rule
for t = 1 : (trials-1)

    betaDist_right(t+1, :) = alphaForget * betaDist_right(t, :) + (1 - alphaForget).*[priorA priorB] ;
    betaDist_left(t+1, :) = alphaForget * betaDist_left(t, :) + (1 - alphaForget).*[priorA priorB] ; % 1 - alphaForget decays it down to priorAlpha, priorBeta
    
    
    if choice(t, 1) == 1 % right choice
        if outcome(t, 1) == 1 % rewarded
            betaDist_right(t+1, 1) = betaDist_right(t+1, 1) + alphaLearn_reward*1; 
        else
            betaDist_right(t+1, 2) = betaDist_right(t+1, 2) + alphaLearn_noreward*1;
        end
    elseif choice(t, 2) == 1 % left choice
        if outcome(t ,2) == 1 % rewarded
            betaDist_left(t+1, 1) = betaDist_left(t+1, 1) + alphaLearn_reward*1;
        else
            betaDist_left(t+1, 2) = betaDist_left(t+1, 2) + alphaLearn_noreward*1;
        end
    end
    
%     uR = betaDist_right(t+1,1) / (betaDist_right(t+1,1)+ betaDist_right(t+1,2));
%     uL = betaDist_left(t+1,1)  / (betaDist_left(t+1,1) + betaDist_left(t+1,2));
%     sigR = (betaDist_right(t+1, 1) * betaDist_right(t+1, 2)) / ((betaDist_right(t+1, 1) + betaDist_right(t+1, 2))^2 * ...
%             (betaDist_right(t+1, 1) + betaDist_right(t+1, 2) + 1));
%     sigL = (betaDist_left(t+1, 1) * betaDist_left(t+1, 2)) / ((betaDist_left(t+1, 1) + betaDist_left(t+1, 2))^2 * ...
%             (betaDist_left(t+1, 1) + betaDist_left(t+1, 2) + 1));
%     Z = ((uR - uL) / ((sigR + sigL)^(1/2)));
%     probChoice(t+1, :) = [normcdf(Z) (1 - normcdf(Z))];

    Ra = betaDist_right(t+1, 1);
    Rb = betaDist_right(t+1, 2);
    La = betaDist_left(t+1, 1);
    Lb = betaDist_left(t+1, 2);
    fun = @(x)(x.^(Ra - 1).*(1 - x).^(Rb - 1))/beta(Ra, Rb).*betacdf(x, La, Lb);
    tP = integral(fun, 0, 1); % probability of right choice

    probChoice(t+1, :) = [tP, (1-tP)];
     
end


% To calculate likelihood:
LH = likelihood(choice, probChoice);
end

