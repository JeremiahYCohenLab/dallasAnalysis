function outputStruct = Qlearning(filename)

%% Determine file path and parse data
% [sessionData, blockSwitch] = loadBehavioralData(fileOrFolder);
% b = parseBehavioralData(sessionData, blockSwitch);

[sessionData, blockSwitch, out] = loadBehavioralData(filename);
b = parseBehavioralData(sessionData, blockSwitch);

%% Q learning model

alpha = 0:.01:1;
beta = 0:.05:20;
likelihood = zeros(length(beta),length(alpha));
for currAlpha = 1:length(alpha)
    for currBeta = 1:length(beta)
        Ql = 0;
        Qr = 0;
        for i = 1:length(b.responseInds)
            if b.allChoices(i) == 1 % right choice
                Qr = Qr + alpha(currAlpha)*(abs(b.allRewards(i)) - Qr);
                likelihood(currBeta,currAlpha) = likelihood(currBeta,currAlpha) + log(exp(beta(currBeta)*Qr)/(exp(beta(currBeta)*Qr)+exp(beta(currBeta)*Ql)));
            elseif b.allChoices(i) == -1 % left choice
                Ql = Ql + alpha(currAlpha)*(abs(b.allRewards(i)) - Ql);
                likelihood(currBeta,currAlpha) = likelihood(currBeta,currAlpha) + log(exp(beta(currBeta)*Ql)/(exp(beta(currBeta)*Qr)+exp(beta(currBeta)*Ql)));
            end
        end
    end
end

[~,idx] = max(likelihood(:)); 
[ml_beta_idx, ml_alpha_idx] = ind2sub(size(likelihood),idx);
ml_alpha = alpha(ml_alpha_idx);
ml_beta = beta(ml_beta_idx);

% Choices
P_rightChoice = [];
Qr = 0;
Ql = 0;
peR = [];
peL = [];
pe = [];
for i = 1:length(b.responseInds)
    P_rightChoice(i) = exp(ml_beta*Qr)/(exp(ml_beta*Qr) + exp(ml_beta*Ql));
    if b.allChoices(i) == 1 % right choice
        Qr = Qr + ml_alpha*(abs(b.allRewards(i)) - Qr);
        peR = [peR (abs(b.allRewards(i)) - Qr)]; 
        pe = [pe (abs(b.allRewards(i)) - Qr)];
    else % left choice
        Ql = Ql + ml_alpha*(abs(b.allRewards(i)) - Ql); 
        peL = [peL (abs(b.allRewards(i)) - Ql)];
        pe = [pe (abs(b.allRewards(i)) - Ql)];
    end
end

%     map = repmat(linspace(25,230,50)./255,3,1)';
%     colormap(map)
%     figure; subplot(2,1,1);
%     imagesc(beta,alpha,likelihood')
%     xlabel('\beta')
%     ylabel('\alpha')
%     set(gca,'YDir','normal')
%     
%     subplot(2,1,2); hold on
%     normKern = normpdf(-15:15,0,2);
%     normKern = normKern / sum(normKern);
%     xVals = (1:(length(normKern) + length(b.allChoices) - 1)) - round(length(normKern)/2);
%     plot(xVals, conv(b.allChoices,normKern)/max(conv(b.allChoices,normKern)),'k','linewidth',2);
%     plot(xVals, conv(b.allRewards,normKern)/max(conv(b.allRewards,normKern)),'--','Color',[100 100 100]./255,'linewidth',2)
%     xlabel('Trials')
%     ylabel('<-- Left       Right -->')
%     legend('Choices','Rewards')
%     xlim([1 length(b.allChoice_R)])
%     ylim([-1 1])
%     for i = 1:length(blockSwitch)
%         bs_loc = blockSwitch(i);
%         plot([bs_loc bs_loc],[-1 1],'r--')
%     end
%     
%     plot(P_rightChoice*2-1,'b--')

outputStruct.likelihood = likelihood;
outputStruct.alpha = ml_alpha;
outputStruct.beta = ml_beta;
outputStruct.peR = peR;
outputStruct.peL = peL;
outputStruct.pe = pe;