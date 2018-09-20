% This function computes a lick density function
% x is a lick train of zeroes and ones (in ms)
% Written by Bilal Bari - 2014-12-10

function myLDF = ldf(x, tc)
    % This is a boxcar filter
    myFilt = ones(1,tc);
    myFilt = 1000 * myFilt./sum(myFilt);
    myLDF = conv(myFilt, x);