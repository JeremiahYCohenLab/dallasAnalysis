% This function computes the area under the receiver operating
%   characteristic curve
% x is a matrix of spike density functions, one for each trial
% binSize is the bin size to use for x
% Analysis is from:
%   Cohen JY, Haesler S, Vong L, Lowell BB, Uchida N. Neuron-type-specific
%     signals for reward and punishment in the ventral tegmental area.
%     Nature 482: 85-88, 2012.
% Written by Jeremiah Y. Cohen

function myAUROC = auROC(x)
    % Bin the spike density function in 100-ms intervals
    binnedSDF = zeros(length(x), 50);
    for i = 1 : length(x)
        temp = [];
        temp2 = x(i, :);
        for j = 1 : 100 : 4901
            temp = cat(1, temp, mean(temp2(j : (j + 99))));
        end
        binnedSDF(i, :) = temp;
    end
    
    % Baseline is the 900 ms before CS onset
    baselineFR = binnedSDF(:, 2 : 9);
    
    % Calculate the ROC curve and the area under the curve
    rocBaseline = [];
    rocMain = [];
    myAUROC = [];
    for i = 1 : 50
        for j = ceil(max(max(binnedSDF))) : -.1 : 0
            rocBaseline = cat(2, rocBaseline, ...
                length(sum(baselineFR(:, i) >= j))/length(baselineFR(:, i)));
            rocMain = cat(2, rocMain, ...
                length(sum(binnedSDF(:, i) >= j))/length(binnedSDF(:, i)));
        end
        myAUROC = cat(1, myAUROC, cumtrapz(rocBaseline(i, :), ...
            rocMain(i, :)));
    end

    