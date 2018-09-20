function compareChoiceByITI_opMD_tmp(xlFile, animals, pre, post, revForFlag)

 if nargin < 5
     revForFlag = 0;
 end
 
 figure;
 l = ceil(length(animals))/2;
 for i = 1:length(animals)

    [x, probSwitchNoRwd_pre, probStayRwd_pre] = combineChoiceByITI_opMD_tmp(xlFile, animals{i}, pre, revForFlag);
    close;
    [x, probSwitchNoRwd_post, probStayRwd_post] = combineChoiceByITI_opMD_tmp(xlFile, animals{i}, post, revForFlag);
    close;
     
    subplot(2,l,i); 
    hold on;
    plot(x, probSwitchNoRwd_pre, '-m', 'linewidth', 2)
    plot(x, probSwitchNoRwd_post, '-', 'Color', [1 0.7 1], 'linewidth', 2)
    plot(x, probStayRwd_pre, '-c', 'linewidth', 2)
    plot(x, probStayRwd_post, '-', 'Color', [0.7 1 1], 'linewidth', 2)
    if i == 1
        legend('P(stay|rwd) pre', 'P(stay|rwd) post', 'P(switch|no rwd) pre', 'P(switch|no rwd) post')
    end
    xlabel('ITI length')
    title(animals{i})
 end
 
  figure;
 l = ceil(length(animals))/2;
 for i = 1:length(animals)
