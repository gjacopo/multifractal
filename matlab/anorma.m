%% ANORMA - Substract the mean of a signal to the signal.
%
%% Syntax
%    [signal,med] = anorma(signal)

function [signal,med] = anorma(signal)

% [sx,sy] = size(signal);
% med = sum(sum(signal)) / (sx*sy);
med = mean(mean(signal))

signal = signal - med;

