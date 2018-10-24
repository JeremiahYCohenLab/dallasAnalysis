function laserInput = desiredLaserInput_intanF(desiredOutput, laserSource, opticFiberTrans)
% laserInput = desiredLaserInput(desiredOutput, opticFiberTrans)
%   DESCRIPTION
%       Uses a calibration range from 1 to 10 fit with output = a*exp(b*input)
%   OUTPUT(S)
%       laserInput: knob setting for a given laser to achive a desired output
%   INPUT(S)
%       desiredOutput: output at most distal optic fiber (mW)
%       laserSource: integer; either '532' or '473'
%       opticFiberTrans: vector of transmission percentage for optic fiber(s)

if nargin < 3
    opticFiberTrans = 0.8;
end

% calibrated on 20180820
if laserSource == 473
    laserInput = ((desiredOutput*2/opticFiberTrans)/ 0.00387)^(1/4.278);
% elseif laserSource == 532
%     p = [0.3323 11.64 (-26.09 - desiredOutput/opticFiberTrans)];
%     r = roots(p);
%     laserInput = max(r);
else
    error('laserSource must be either 473 or 532')
end