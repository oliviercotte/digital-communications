% File: sim4PAM_AWGN_SQRC.m
% Simulation of 4-PAM mapping at baseband over AWGN channel using a 
% Square-Root-Raised-Cosine filter as pulse shaping filter and
% matched filter.
 
%% Close everything, reset workspace
close all; % Closes all of the figures that you have generated in your program
clear all; % Deletes all stored variables in your workspace 
clc; % Removes all lines in your command window

%% Creates a different seed each time
rng shuffle;
 
%% Simulation parameters
N = 2 ^ 18; % Number of symbols
M = 4; % Size of signal constellation
k = log2(M); % Number of bits per symbol
NumSamplesPerSymbol = 4; % Oversampling factor
 
%% Square-Root-Raised-Cosine filter filter parameters
Span = 10; % Filter span in symbols
Rolloff = 0.25; % Rolloff factor of filter
Syms = -Span : 1 / NumSamplesPerSymbol : Span;
NumCoefficients = (2 * Span * NumSamplesPerSymbol) + 1; % Number of coefficients
PSF = sqrc(Syms, Rolloff); % Square-Root-Raised-Cosine filter
 
%% Normalize psf such that its energy is 1
% To make sure the recevied signal energy does not depend on the psf
EnergyPSF = sum(PSF .^ 2);
PSF = PSF ./ sqrt(EnergyPSF);
 
%% PAM4 mapping (Generation of equiprobable input symbols)
% Generate random 4-PAM (+/-1, +/-3) sequence
% alpha4pam = [-3 -1 1 3]; % 4-PAM alphabets
% Message = randsrc(1, N, alpha4pam);
Message = m_pam(N, M);

%% Normalization of signal energy to 1
TransmittedsSignal = Message ./ sqrt(5);
 
%% SNR (Es/No) values
Eb_No_dB = -10 : 1 : 10; % Signal to Noise Ratio in dB
Eb_No = 10 .^ (Eb_No_dB ./ 10); % Signal to Noise Ratio in Linear 
Es_No_dB  = Eb_No_dB + 10 .* log10(k);
Es_No = 10 .^ (Es_No_dB ./ 10);
 
%% Generation of Noise (mean = 0, std = sqrt(No / 2)
No = 10 .^ (-Es_No_dB ./ 10); % AWGN single-sided power
Sigma = sqrt(No ./ 2); % AWGN standard deviation
 
%% Upsampling
UpsampledData = upsample(TransmittedsSignal, NumSamplesPerSymbol);
 
%% Pulse shaping filter
PSFSignal = conv(PSF, UpsampledData);
 
%% AWGN CHANNEL
 % additive white gaussian noise
AWGN = randn(1, length(PSFSignal));
Noise = Sigma' * AWGN;
DataSent = cellfun(@(x) x + PSFSignal, num2cell(Noise, 2),...
    'UniformOutput', false);
 
%% Matched Filter
MatchedFilter = cellfun(@(tx) conv(PSF, tx), DataSent,...
    'UniformOutput', false);
 
%% Downsampling
FullLength = length(PSF) + length(AWGN) - 1;
MatchedFilterOutput = cell2mat(cellfun(@(x) x(NumCoefficients : NumSamplesPerSymbol : FullLength - NumCoefficients),...
    MatchedFilter, 'UniformOutput', false));
 
%% PAM4 constellation
lambda1 = -2/sqrt(5); lambda2 = 0; lambda3 = 2/sqrt(5);
 
%% DECISION DEVICE FOR 4-PAM
EstimatedMessage(MatchedFilterOutput < lambda1) = -3;
EstimatedMessage(MatchedFilterOutput >= lambda1 & MatchedFilterOutput < lambda2) = -1;
EstimatedMessage(MatchedFilterOutput >= lambda2 & MatchedFilterOutput < lambda3) = 1;
EstimatedMessage(MatchedFilterOutput >= lambda3) = 3;
EstimatedMessage = reshape(EstimatedMessage, [], N);
 
%% COMPUTE BIT ERRORS FOR 4-PAM
SimulationBER = cell2mat(cellfun(@(x) size(find(Message - x), 2),...
    num2cell(EstimatedMessage, 2), 'UniformOutput', false))' ./ (k * N);

%% SHOW THE PLOT
% After the simulation over different SNR values, a vector of BER is
% obtained with  respect  to  the  SNR  vector  previously  defined.
figure

%% Plot theoretical value of bit error probability for 4-PAM (with SRRC)
TheoriticalBERWithSQRC = (3 / 4) * qfunc(sqrt((2 / 5) * Es_No));
semilogy(Eb_No_dB, TheoriticalBERWithSQRC, 'bs','LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',8), hold on;

%% %% Plot simulation BER
semilogy(Eb_No_dB, SimulationBER, 'b--', 'LineWidth', 1.50);

grid on, xlabel('E_b/N_0 (dB)'), ylabel('Bit Error Rate')
legend('Theory 4-PAM with SRRC', 'Simulation 4-PAM with SRRC')
title('Simulation of 4-PAM over AWGN channel using a Square-Root-Raised-Cosine filter as pulse shaping filter'), hold off

%% Print theoritical BER vs simulated
fprintf('\nSimulation of 4-PAM at baseband over AWGN channel\n')
fprintf('using a Square-Root-Raised-Cosine filter as pulse shaping\n') 
fprintf('filter and matched filter.\n');
fprintf('--------------------------------------\n');
fprintf('Eb/No\t\tSim(SRRC)\tTheo(SRRC)')
fprintf('\n--------------------------------------\n');
arrayfun(@(w,x,y) fprintf('%5.2f \t %e\t %e\t\n',w,x,y), Eb_No_dB,...
    SimulationBER, TheoriticalBERWithSQRC);
fprintf('--------------------------------------\n');