close all
clear
clc

addpath functions

%% Parameter setup
% Iteration parameters
SP.iter             = 10;
SP.SNR_db_array     = 0; %-10:5:15;
% OMP parameters
SP.algorithms       = {'OMP'}; % {'LS', 'OMP', 'BPDN', 'IHT'};
SP.criteria         = 'residual'; % 'paths' or 'residual'
% System parameter[s
SP.Nt               = 64;
SP.Nr               = 16;
SP.Lt               = 8;
SP.sym              = 8; % Mx
SP.Lr               = 4;
SP.gridMultiplier   = 1.5;
SP.Gt               = SP.Nt * SP.gridMultiplier;
SP.Gr               = SP.Nr * SP.gridMultiplier;
SP.Mmax             = SP.Nt*SP.Nr/SP.Lr;
% SP.permIter         = 1;
SP.bandwidth        = 500e6;
SP.sigma            = sqrt(10^(.1*(-174+10*log10(SP.bandwidth))));
SP.rfArchitecture   = 'PS'; 
% Following two parameters only applicable when 'PS' architecture is used.
SP.rfBFtype         = 'DFT'; % RF beamformer type. 'DFT', 'DCT', 'Hadamard', 'GD'
SP.rfQuant          = true; % RF precoder/combiner quantization (HW constraint).
SP.rfQuantBits      = 6;
% Channel-related parameters
SP.Np               = 4; % # of total paths

%% Simulation
SP = main_function(SP);
