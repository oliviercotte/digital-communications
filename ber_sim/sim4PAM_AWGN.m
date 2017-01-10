% File: sim4PAM_AWGN.m
% Simulation of 4-PAM mapping at baseband over AWGN channel
% Sketch theoretical and simulated bit error rate of the system versus
% Eb/No.

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

%% PAM4 mapping (Generation of equiprobable input symbols)
% Generate random 4-PAM (+/-1, +/-3) sequence
% alpha4pam = [-3 -1 1 3]; % 4-PAM alphabets
% Message = randsrc(1, N, alpha4pam);
Message = m_pam(N, M);

%% Normalization of signal energy to 1
TransmittedsSignal = Message ./ sqrt(5);

%% SNR (Es / No) values
Eb_No_dB = -10 : 1 : 10; % Signal to Noise Ratio in dB
Eb_No = 10 .^ (Eb_No_dB ./ 10); % Signal to Noise Ratio in Linear 
Es_No_dB  = Eb_No_dB + 10 .* log10(k);
Es_No = 10 .^ (Es_No_dB ./ 10);

%% Generation of Noise N (mean = 0, std = sqrt(No / 2)
No = 10 .^ (-Es_No_dB ./ 10); % AWGN single-sided power
Sigma = sqrt(No ./ 2); % AWGN standard deviation

%% AWGN CHANNEL
AWGN = randn(1, N); % additive white gaussian noise
Noise = Sigma' * AWGN;
ReceivedSignal = cell2mat(cellfun(@(x) x + TransmittedsSignal,...
    num2cell(Noise, 2), 'UniformOutput', false));

%% PAM4 constellation
lambda1 = -2/sqrt(5); lambda2 = 0; lambda3 = 2/sqrt(5);

%% DECISION DEVICE FOR 4-PAM
EstimatedMessage(ReceivedSignal < lambda1) = -3;
EstimatedMessage(ReceivedSignal >= lambda1 & ReceivedSignal < lambda2) = -1;
EstimatedMessage(ReceivedSignal >= lambda2 & ReceivedSignal < lambda3) = 1;
EstimatedMessage(ReceivedSignal >= lambda3) = 3;
EstimatedMessage = reshape(EstimatedMessage, [], N);
    
%% COMPUTE BIT ERRORS FOR 4-PAM
SimulationBER = cell2mat(cellfun(@(x) size(find(Message - x), 2),...
    num2cell(EstimatedMessage, 2), 'UniformOutput', false))' ./ (k * N);

%% SHOW THE PLOT
% After the simulation over different SNR values, a vector of BER is
% obtained with  respect  to  the  SNR  vector  previously  defined.
figure

%% Plot theoretical value of bit error probability for 4-PAM
TheoriticalBER = (3 / 4) * qfunc(sqrt((2 / 5) * Es_No));
semilogy(Eb_No_dB, TheoriticalBER, 'bs','LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',8), hold on;

%% Plot simulation BER
semilogy(Eb_No_dB, SimulationBER, 'b--', 'LineWidth', 1.50);

grid on, xlabel('E_b/N_0 (dB)'), ylabel('Bit Error Rate')
legend('Theory 4-PAM', 'Simulation 4-PAM')
title('Simulation of 4-PAM at baseband over AWGN channel')
hold off

%% Print theoritical BER vs simulated
fprintf('\nSimulation of 4-PAM at baseband over AWGN channel\n');
fprintf('Eb/No \tSIM_4PAM \t\tTHEO_4PAM\n--------------------------------------\n');
arrayfun(@(x,y,z) fprintf('%5.2f \t%e\t%e\n',x,y,z), Eb_No_dB,...
    SimulationBER, TheoriticalBER);
fprintf('--------------------------------------\n');