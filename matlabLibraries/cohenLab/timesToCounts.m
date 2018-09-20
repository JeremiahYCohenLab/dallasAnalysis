% This function changes from spike times to counts
% Right now, this is only implemented for trials of length 5001 and bin
%   size = 1 ms, for data within a 3D array
% Written by Jeremiah Y. Cohen

function myData = timesToCounts(x)
    myData = zeros(1, 5001);
    temp = x(1, 1, :);
    myData(temp(~isnan(temp)) + 1000) = 1;
    