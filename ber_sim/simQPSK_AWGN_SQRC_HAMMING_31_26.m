% File: simQPSK_AWGN_SQRC_HAMMING_31_26
% Apply Hamming (31, 26) as channel coding to the QPSK system.

%% Close everything, reset workspace
clear all;
close all;
clc;

%% Creates a different seed each time
rng shuffle;

%% Generator Matrix G, and Parity Matrix H
% We want a [ 31, 26 ] Hamming code.
m = 5;
n_code = 2 ^ m - 1; % % Codeword (Block) length
k_code = n_code - m; % Message size
t_code = hammingbound(n_code, k_code); % Error Correction Capability
Rc = k_code / n_code ;
%[H,G] = hammgen(m); % Produce the parity-check matrix.
P = [0 0 0 1 1;0 0 1 0 1;0 0 1 1 0;0 1 0 0 1;0 1 0 1 0;0 1 1 0 0;1 0 0 0 1;1 0 0 1 0;1 0 1 0 0;1 1 0 0 0;0 0 1 1 1;0 1 0 1 1;0 1 1 0 1;0 1 1 1 0;1 0 0 1 1;1 0 1 0 1;1 0 1 1 0;1 1 0 0 1;1 1 0 1 0;1 1 1 0 0;0 1 1 1 1;1 0 1 1 1;1 1 0 1 1;1 1 1 0 1;1 1 1 1 0;1 1 1 1 1];
G = [P, eye(k_code)];
H = [eye(n_code - k_code), P'];ErrorPatterns = combis(n_code, t_code); % Error patterns (Fig.9.8)
SyndromeMatrix = rem(ErrorPatterns * H', 2); % The syndrome matrix

%% Simulation parameters
N = 2 ^ 15; % Number of codewords
M = 4; % Conselation size
k = log2(M); % Bits per symbol
NumSamplesPerSymbol = 4; % Oversampling factor
BitLength = N * k_code / Rc;

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

%% SNR (Es / No) values
Eb_No_dB = -10 : 1 : 8; % Signal to Noise Ratio in dB
Eb_No = 10 .^ (Eb_No_dB ./ 10);
Es_No_dB = Eb_No_dB + 10 .* log10(k);
Es_No = 10 .^ (Es_No_dB ./ 10);
Ec_No_dB = Eb_No_dB + 10 .* log10(k * Rc);
Ec_No = 10 .^ (Ec_No_dB ./ 10);
sqrtEc_No = sqrt(Ec_No);

%% Create random input message block (N*k)
% Message vector
Message_I = floor(rand(N, k_code) .* 2);
Message_Q = floor(rand(N, k_code) .* 2);

%% Encoding c = mG (mod 2)
% Coded vector
Codeword_I = cell2mat(cellfun(@(x) rem(x * G, 2), num2cell(Message_I, 2), 'UniformOutput', false));
Codeword_Q = cell2mat(cellfun(@(x) rem(x * G, 2), num2cell(Message_Q, 2), 'UniformOutput', false));

%% Simulate QPSK modulator
InPhase = (2 * Codeword_I - 1); %In-phase symbol generation
Quadrature = (2 .* Codeword_Q - 1); %Quadrature symbol 

%% Upsampling
data_I_sampled = cellfun(@(x) upsample(x, NumSamplesPerSymbol),...
    num2cell(InPhase, 2), 'UniformOutput', false);
data_Q_sampled = cellfun(@(x) upsample(x, NumSamplesPerSymbol),...
    num2cell(Quadrature, 2), 'UniformOutput', false);

%% Pulse shaping filter
data_I_sent = cell2mat(cellfun(@(x) conv(PSF, x), data_I_sampled,...
    'UniformOutput',false));
data_Q_sent = cell2mat(cellfun(@(x) conv(PSF, x), data_Q_sampled,...
    'UniformOutput', false));

%% Generation of Noise (mean = 0, std = sqrt(No / 2)
% Additive white gaussian noise
LengthNoise = NumSamplesPerSymbol * n_code + length(PSF) - 1;
AWGN_I = randn(N, LengthNoise);
AWGN_Q = randn(N, LengthNoise);
Noise_I = cellfun(@(x) AWGN_I ./ x, num2cell(sqrtEc_No), 'UniformOutput', false);
Noise_Q = cellfun(@(x) AWGN_Q ./ x, num2cell(sqrtEc_No), 'UniformOutput', false);

%% AWGN CHANNEL
% Received vector with noise
data_I_received = cell2mat(cellfun(@(x) data_I_sent + x, Noise_I, 'UniformOutput', false));
data_Q_received = cell2mat(cellfun(@(x) data_Q_sent + x, Noise_Q, 'UniformOutput', false));
data_I_received = reshape(data_I_received, N, LengthNoise, length(Eb_No_dB));
data_Q_received = reshape(data_Q_received, N, LengthNoise, length(Eb_No_dB));

%% Matched Filter
data_I_deMod = cellfun(@(tx) conv(PSF, tx), num2cell(data_I_received, 2),...
    'UniformOutput', false);
data_Q_deMod = cellfun(@(tx) conv(PSF, tx), num2cell(data_Q_received, 2),...
    'UniformOutput', false);

%% Downsampling
FullLength = length(PSF) + LengthNoise - 1;
data_I_deMod = cell2mat(cellfun(@(x) x(NumCoefficients : NumSamplesPerSymbol : FullLength - NumCoefficients),...
    data_I_deMod, 'UniformOutput', false));
data_Q_deMod = cell2mat(cellfun(@(x) x(NumCoefficients : NumSamplesPerSymbol : FullLength - NumCoefficients),...
    data_Q_deMod, 'UniformOutput', false));

%% DEMODULATION
% Hard quantize demodulator
% All positive received bits converted to +1
% All negative received bits converted to 0
data_I_sliced = data_I_deMod > 0; % In-phase demodulation
data_Q_sliced = data_Q_deMod > 0; % Quadrature demodulation

%% Compute syndrome
Syndrome_I = cellfun(@(x) rem(x * H',2), num2cell(data_I_sliced, [1 2]), 'UniformOutput', false);
Syndrome_Q = cellfun(@(x) rem(x * H',2), num2cell(data_Q_sliced, [1 2]), 'UniformOutput', false);

%% Find estimated error pattern in syndrome table
Correction_I = cellfun(@(x) ComputeEstimatedErrorPatternFromSyndrome(x, SyndromeMatrix, ErrorPatterns), Syndrome_I, 'UniformOutput', false);
Correction_Q = cellfun(@(x) ComputeEstimatedErrorPatternFromSyndrome(x, SyndromeMatrix, ErrorPatterns), Syndrome_Q, 'UniformOutput', false);

%% Compute corrected codeword
ReceivedCodeword_I = cellfun(@(x, y) rem(x + y, 2), num2cell(data_I_sliced, [1 2]), Correction_I, 'UniformOutput', false);
ReceivedCodeword_Q = cellfun(@(x, y) rem(x + y, 2), num2cell(data_Q_sliced, [1 2]), Correction_Q, 'UniformOutput', false);

%% Extract the estimated message
EstimatedMessage_I = cellfun(@(x) x(:, n_code - k_code + 1: n_code), ReceivedCodeword_I, 'UniformOutput', false);
EstimatedMessage_Q = cellfun(@(x) x(:, n_code - k_code + 1: n_code), ReceivedCodeword_Q, 'UniformOutput', false);

%% Computes simulated bit error rate
% Message bit error probability;
Error_I = cell2mat(cellfun(@(x) sum(sum(Message_I ~= x)), EstimatedMessage_I, 'UniformOutput', false)); %In-phase BER calculation
Error_Q = cell2mat(cellfun(@(x) sum(sum(Message_Q ~= x)), EstimatedMessage_Q, 'UniformOutput', false)); %Quadrature BER calculation
Error = reshape(mean([Error_I Error_Q]), 1, length(Eb_No_dB));
BER = Error ./ BitLength;  %Overall BER

%% Theoritical Bit Error Prob. for QPSK w and w/o coding
Pb_Uncoded = qfunc(sqrt(Es_No)); % Uncoded msg BER with BPSK by Eq.(7.3.4)
Pc = qfunc(sqrt(Ec_No)); % Crossover probability by Eq.(7.3.4)
Pb = prob_err_msg_bit(Pc, n_code, t_code);

%% SHOW THE PLOT
% After the simulation over different SNR values, a vector of BER is obtained
% with  respect  to  the  SNR  vector  previously  defined.
figure

%% Plot theoretical value of bit error probability for QPSK
semilogy(Eb_No_dB, Pb,  'bs','LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',8), hold on;

%% Plot simulation BER of the modem
semilogy(Eb_No_dB, Pc, 'r-','LineWidth', 2.0)

%% Plot simulation BER of QPSK without coding
semilogy(Eb_No_dB, Pb_Uncoded, 'm-<', 'linewidth', 2.0);

%% Plot simulation BER of QPSK with coding
semilogy(Eb_No_dB, BER, 'b-','LineWidth', 2.0)

%% Title, legend, etc
grid on, axis tight, xlabel('E_b/N_0 (dB)'), ylabel('Bit Error Rate')
title('Simulation of QPSK modulation scheme at baseband over AWGN channel with and without coding')
legend('Theory BER of QPSK with Hamming (31, 26)', 'P_c of the modem', 'Theory uncoded QPSK', ...
    'Simulated BER of QPSK with Hamming (31, 26)'), grid on, hold off;