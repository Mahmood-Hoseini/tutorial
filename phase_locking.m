% copyright (c) 2017 Mahmood Hoseini

% This program is free software: you can redistribute it and/or modify it
% under the terms of the gnu general public license as published by the
% free software foundation, either version 3 of the license, or (at your
% option) any later version. you should have received a copy of the gnu
% general public license along with this program.  if not,
% see <http://www.gnu.org/licenses/>.

clc; clear all;
set(0, 'defaultaxesfontsize', 15)
set(0, 'defaulttextfontsize', 12)
set(0, 'defaultaxeslinewidth', 2)


load('/Users/idl/Downloads/KH63_for_Mahmood/KH63_cells.mat')
load('/Users/idl/Downloads/KH63_for_Mahmood/block37_gamma.mat')

freq_bands = [1.0, 2.2;
              2.3, 4.3;
              4.8, 9.3;
              10.0, 24;
              25.0, 49;
              52.0, 100];


%% Parameters
cell_id = 24;
chan = 11;
Fs = 400; % Sampling rate in Hz
dt = 1/Fs; % Time step size
pad = 1;      % pad the time series with zeroes (recommended)
dj = 0.25;    % this will do 4 sub-octaves per octave
s0 = 4*dt;    % this says start at a scale of 6 months
j1 = 7/dj;    % this says do 7 powers-of-two with dj sub-octaves each
lag = 0.72;  % lag-1 autocorrelation for red noise background
mother = 'Morlet';

LFP = LFPdat(chan).data.signal;
n = length(LFP);
time = (1 : n)*dt;  % construct time array

time_diff = unixtime(LFPdat(chan).data.starttime) - unixtime(CELL(chan).EXPTSTART);
relativeStartTime = (((time_diff(3)*24 + time_diff(4))*60 + time_diff(5))*60 + time_diff(6) + CELL(chan).trem);

%% Z-scoreing!
variance = std(LFP)^2;
LFP = (LFP - mean(LFP))/sqrt(variance) ;

%% Wavelet transform:
[wave,  period, scale, coi] = wavelet(LFP, dt, pad, dj, s0, j1, mother);

%% Phase consistency (PC) index
phase_diff = angle(wave);
R_ = cos(phase_diff);
r = zeros(1, length(freq_bands), 2);
f = 1./period1;
for ii = 1 : length(freq_bands)
    choice = find(f > freq_bands(ii, 1) & f < freq_bands(ii, 2));
    r(1, ii, 1) = mean(mean(R_(choice, 1:2*Fs)));
    r(1, ii, 2) = mean(mean(R_(choice, 2*Fs:end)));
end
R_lst = [R_lst; r];