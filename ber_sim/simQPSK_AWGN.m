% File: simQPSK_AWGN.m
% Simulation of Quadrature-Phase-Shift-Keying (QPSK) modulation scheme at
% baseband over AWGN channel. Sketch theoretical and simulated bit error 
% rate of the system versus Eb/No.

%% Close everything, reset workspace
close all; % Closes all of the figures that you have generated in your program
clear all; % Deletes all stored variables in your workspace 
clc; % Removes all lines in your command window

%% Creates a different seed each time
rng shuffle;

%% Simulation parameters
N = 2 ^ 20; % Number of symbols
M = 4; % Size of signal constellation
k = log2(M); % Number of bits per symbol

%% QPSK mapping (Generation of equiprobable input symbols)
% Generate random QPSK (+/-1, +/-1) sequence
Message_I = 2 * (round(rand(1, N)) - 0.5);
Message_Q = 2 * (round(rand(1, N)) - 0.5);

%% Normalization of signal energy to 1
InPhase = Message_I ./ sqrt(2);
Quadrature = Message_Q ./ sqrt(2);

%% SNR (Es / No) values
Eb_No_dB = -10 : 1 : 10; % Signal to Noise Ratio in dB
Eb_No = 10 .^ (Eb_No_dB ./ 10); % Signal to Noise Ratio in Linear 
Es_No_dB  = Eb_No_dB + 10 .* log10(k);
Es_No = 10 .^ (Es_No_dB ./ 10);

%% Generation of Noise N (mean = 0, std = sqrt(No / 2)
No = 10 .^ (-Es_No_dB ./ 10); % AWGN single-sided power
Sigma = sqrt(No ./ 2); % AWGN standard deviation

%% AWGN CHANNEL
% Additive white gaussian noise
AWGN_I = randn(1, N);
AWGN_Q = randn(1, N);
Noise_I = Sigma' * AWGN_I;
Noise_Q = Sigma' * AWGN_Q;
InPhaseReceiver  = cell2mat(cellfun(@(x) x + InPhase,...
    num2cell(Noise_I, 2), 'UniformOutput', false));
QuadratureReceiver  = cell2mat(cellfun(@(x) x + Quadrature,...
    num2cell(Noise_Q, 2), 'UniformOutput', false));

%% DEMODULATION
EstimatedMessage_I = sign(InPhaseReceiver);
EstimatedMessage_Q = sign(QuadratureReceiver);

%% COMPUTE BIT ERRORS FOR QPSK
SimulationBER_I = cellfun(@(x) N - sum(Message_I == x),...
    num2cell(EstimatedMessage_I, 2), 'UniformOutput', false); 
SimulationBER_Q = cellfun(@(x) N - sum(Message_Q == x),...
    num2cell(EstimatedMessage_Q, 2), 'UniformOutput', false);
SimulationBER = cell2mat(cellfun(@(x,y) mean([x y]),...
    SimulationBER_I, SimulationBER_Q,...
    'UniformOutput', false))' ./ N;

%% SHOW THE PLOT
% After the simulation over different SNR values, a vector of BER is
% obtained with  respect  to  the  SNR  vector  previously  defined.
figure

% Plot theoretical value of bit error probability for QPSK
TheoriticalBER = qfunc(sqrt(Es_No));
semilogy(Eb_No_dB, TheoriticalBER, 'bs', 'LineWidth', 2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',8), hold on;

% Plot simulation BER of QPSK
semilogy(Eb_No_dB, SimulationBER, 'b--','LineWidth', 1.50)

grid on, xlabel('E_b/N_0 (dB)'), ylabel('Bit Error Rate')
title('Simulation of QPSK with transmission over an AWGN channel')
legend('Theory QPSK', 'Simulation QPSK')
grid on
hold off

fprintf('\nSimulation of QPSK with transmission over an AWGN channel\n');
fprintf('--------------------------------------\n');
fprintf('Eb/No  \tBER_4PAM\n--------------------------------------\n');
fprintf('Eb/No \tSIM_QPSK \t\tTHEO_QPSK\n--------------------------------------\n');
arrayfun(@(x,y,z) fprintf('%5.2f \t%e\t%e\n',x,y,z), Eb_No_dB,...
    SimulationBER, TheoriticalBER);
fprintf('--------------------------------------\n');