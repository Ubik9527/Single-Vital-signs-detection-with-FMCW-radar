function Hd = IIR_breath
%IIR_BREATH 返回离散时间滤波器对象。

% MATLAB Code
% Generated by MATLAB(R) 9.13 and Signal Processing Toolbox 9.1.
% Generated on: 20-Nov-2023 14:30:34

% Butterworth Bandpass filter designed using FDESIGN.BANDPASS.

% All frequency values are in Hz.
Fs = 50;  % Sampling Frequency

N   = 4;    % Order
Fc1 = 0.15;  % First Cutoff Frequency
Fc2 = 0.4;  % Second Cutoff Frequency

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, Fs);
Hd = design(h, 'butter');

% [EOF]
