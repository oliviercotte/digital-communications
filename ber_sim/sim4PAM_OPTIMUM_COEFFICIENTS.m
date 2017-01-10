% File: sim4PAM_OPTIMUM_COEFFICIENTS.m
% This script the optimum number of coefficients required per each filter 
% for 4PAM mapping such that the overall bit error rate is close to theory.

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

%% SHOW THE PLOT
% After the simulation over different SNR values, a vector of BER is
% obtained with  respect  to  the  SNR  vector  previously  defined.
figure

TheoriticalBERWithSQRC = (3 / 4) * qfunc(sqrt((2 / 5) * Es_No));

fprintf('\nSimulation of 4-PAM at baseband over AWGN channel\n')
fprintf('using a Square-Root-Raised-Cosine filter as pulse shaping\n')
fprintf('filter and matched filter.\n');
fprintf('--------------------------------------\n');

for ii = 1 : 7
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
    
    %% Euclidean distances
    Dist(ii,:) = pdist2(SimulationBER, TheoriticalBERWithSQRC);
    
    %% Print theoritical BER vs simulated
    fprintf('Eb/No Sim(Srrc span = %d) Theo(Srrc span = %d)', Span, Span)
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

%% Plot theoretical value of bit error probability for 4-PAM (with SRRC)
semilogy(Eb_No_dB, TheoriticalBERWithSQRC, 'bs', 'LineWidth', 1.25, ...
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
title('Simulation of 4-PAM over AWGN channel using a Square-Root-Raised-Cosine filter as pulse shaping filter'), hold off