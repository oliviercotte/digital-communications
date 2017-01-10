% File: simQPSK_OPTIMUM_COEFFICIENTS.m
% This script the optimum number of coefficients required per each filter 
% for QPSK mapping such that the overall bit error rate is close to theory.

%% Close everything, reset workspace
%close all; % Closes all of the figures that you have generated in your program
clear all; % Deletes all stored variables in your workspace
clc; % Removes all lines in your command window

%% Creates a different seed each time
rng shuffle;

%% Simulation parameters
N = 2 ^ 20; % Number of symbols
M = 4; % Size of signal constellation
k = log2(M); % Number of bits per symbol
NumSamplesPerSymbol = 4; % Oversampling factor

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

%% Upsampling
data_I_sampled = upsample(InPhase, NumSamplesPerSymbol);
data_Q_sampled = upsample(Quadrature, NumSamplesPerSymbol);

%% SHOW THE PLOT
% After the simulation over different SNR values, a vector of BER is
% obtained with  respect  to  the  SNR  vector  previously  defined.
figure

TheoriticalBERWithSQRC = qfunc(sqrt(Es_No));

fprintf('\nSimulation of QPSK at baseband over AWGN channel\n')
fprintf('using a Square-Root-Raised-Cosine filter as pulse shaping\n')
fprintf('filter and matched filter.\n');
fprintf('--------------------------------------\n');

for ii = 1 : 10
    %% Square-Root-Raised-Cosine filter filter parameters
    Span = ii; % Filter span in symbols
    Rolloff = 0.25; % Rolloff factor of filter
    Syms = -Span : 1 / NumSamplesPerSymbol : Span;
    NumCoefficients = (2 * Span * NumSamplesPerSymbol) + 1; % Number of coefficients
    PSF = srrc(Syms, Rolloff); % Square-Root-Raised-Cosine filter
    
    %% Normalize psf such that its energy is 1
    % To make sure the recevied signal energy does not depend on the psf
    EnergyPSF = sum(PSF .^ 2);
    PSF = PSF ./ sqrt(EnergyPSF);
    
    %% Pulse shaping filter
    data_I_sent = conv(PSF, data_I_sampled);
    data_Q_sent = conv(PSF, data_Q_sampled);
    
    %% AWGN CHANNEL
    AWGN_I = randn(1, length(data_I_sent));
    AWGN_Q = randn(1, length(data_Q_sent));
    Noise_I = Sigma' * AWGN_I;
    Noise_Q = Sigma' * AWGN_Q;
    data_I_received = cellfun(@(x) x + data_I_sent, num2cell(Noise_I, 2),...
        'UniformOutput', false);
    data_Q_received = cellfun(@(x) x + data_Q_sent, num2cell(Noise_Q, 2),...
        'UniformOutput', false);
    
    %% Matched Filter
    data_I_deMod = cellfun(@(tx) conv(PSF, tx), data_I_received,...
        'UniformOutput', false);
    data_Q_deMod = cellfun(@(tx) conv(PSF, tx), data_Q_received,...
        'UniformOutput', false);
    
    %% Downsampling
    FullLength = length(PSF) + length(AWGN_I) - 1;
    data_I_deMod = cell2mat(cellfun(@(x) x(NumCoefficients : NumSamplesPerSymbol : FullLength - NumCoefficients),...
        data_I_deMod, 'UniformOutput', false));
    data_Q_deMod = cell2mat(cellfun(@(x) x(NumCoefficients : NumSamplesPerSymbol : FullLength - NumCoefficients),...
        data_Q_deMod, 'UniformOutput', false));
    
    %% DEMODULATION
    EstimatedMessage_I = sign(data_I_deMod);
    EstimatedMessage_Q = sign(data_Q_deMod);
    
    %% COMPUTE BIT ERRORS FOR QPSK
    SimulationBER_I = cellfun(@(x) N - sum(Message_I == x),...
        num2cell(EstimatedMessage_I, 2), 'UniformOutput', false);
    SimulationBER_Q = cellfun(@(x) N - sum(Message_Q == x),...
        num2cell(EstimatedMessage_Q, 2), 'UniformOutput', false);
    SimulationBER = cell2mat(cellfun(@(x,y) mean([x y]),...
        SimulationBER_I, SimulationBER_Q,...
        'UniformOutput', false))' ./ N;
    
    %% Euclidean distances
    Dist(ii,:) = pdist2(SimulationBER, TheoriticalBERWithSQRC);
    
    %% Print theoritical BER vs simulated
    fprintf('Eb/No Sim(Srrc span = %d) Theo(SQR-RC width = %d)', Span, Span)
    fprintf('\n--------------------------------------\n');
    arrayfun(@(w,x,y) fprintf('%5.2f \t %e\t %e\t\n',w,x,y), Eb_No_dB,...
        SimulationBER, TheoriticalBERWithSQRC);
    fprintf('--------------------------------------\n');
    
    if Span == 1
        semilogy(Eb_No_dB, SimulationBER, ':p', 'LineWidth', 1.25);
        hold on;
    elseif Span == 2
        semilogy(Eb_No_dB, SimulationBER, ':o', 'LineWidth', 1.25);
    elseif Span == 3
        semilogy(Eb_No_dB, SimulationBER, '--rd', 'LineWidth', 1.25);
    elseif Span == 4
        semilogy(Eb_No_dB, SimulationBER, '--rs', 'LineWidth', 1.25);
    elseif Span == 5
        semilogy(Eb_No_dB, SimulationBER, '--rv', 'LineWidth', 1.25);
    elseif Span == 6
        semilogy(Eb_No_dB, SimulationBER, '-.x', 'LineWidth', 1.25);
    else
        semilogy(Eb_No_dB, SimulationBER, '--rs', 'LineWidth', 1.25);
    end
end

%% Plot theoretical value of bit error probability for QPSK (with SRRC)
semilogy(Eb_No_dB, TheoriticalBERWithSQRC, 'bs','LineWidth', 1.25, ...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',8), hold on;

if Span == 1
    legend('SQR-RC width = 1', 'Theory');
elseif Span == 2
    legend('SQR-RC width = 1', 'SQR-RC width = 2', 'Theory');
elseif Span == 3
    legend('SQR-RC width = 1', 'SQR-RC width = 2', 'Srrc span = 3', 'Theory');
elseif Span == 4
    legend('SQR-RC width = 1', 'SQR-RC width = 2', 'SQR-RC width = 3', 'SQR-RC width = 4', 'Theory');
elseif Span == 5
    legend('SQR-RC width = 1', 'SQR-RC width = 2', 'SQR-RC width = 3', 'SQR-RC width = 4', 'SQR-RC width = 5', 'Theory');
elseif Span == 6
    legend('SQR-RC width = 1', 'SQR-RC width = 2', 'SQR-RC width = 3', 'SQR-RC width = 4', 'SQR-RC width = 5', 'SQR-RC width = 6', 'Theory');
else %Span > 6
    legend('SQR-RC width = 1', 'SQR-RC width = 2', 'SQR-RC width = 3', 'SQR-RC width = 4', 'SQR-RC width = 5', 'SQR-RC width = 6', 'SQR-RC width > 6', 'Theory');
end

grid on, xlabel('E_b/N_0 (dB)'), ylabel('Bit Error Rate')
title('Simulation of QPSK at baseband over AWGN channel using a Square-Root-Raised-Cosine filter as pulse shaping'), hold off