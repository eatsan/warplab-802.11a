
function [output] = channels(input,channel,SNR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% Channels.m
%
% run OFDM time domain samples through channels
%
% Input: %
% input: OFDM time domain samples
% channel: either awgn or Rayleigh as text string
% SNR: signal to noise ratio
% Output: %
% output: channel filtered time domain samples
%
% Author: Keith Howland
% Created: 03 May 07
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%=================================%
%send signal through channel
%=================================%
if strcmp('awgn',channel)
    %Additive White Gaussian Noise Channel
    output = awgn(input,SNR,'measured');
elseif strcmp('rayleigh',channel)
    %Rayleigh Channel with Additive White Gaussian Noise
    chan = rayleighchan(1/(20*10^6),10);
    chan.pathdelays=[0 5e-8];
    chan.avgpathgaindb=[0 0];
    y = filter(chan,input);
    output=awgn(y,SNR,'measured');
elseif strcmp('802_11n_B',channel)
   
    % IEEE 802.11n Channel Model B for SIMO
    
    % Rs: Input Sample Rate
    Rs = 20e6;
    % Nt: Number of transmit antennas
    Nt=2;
    % Nr: Number of receive antennas
    Nr=2;

    tau = [0 10 20 30 40 50 60 70 80] * 1e-9; % Path delays, in seconds
    
    % Average path gains of cluster 1, in dB
    pdb1 = [0 -5.4 -10.8 -16.2 -21.7 -inf -inf -inf -inf];
    % Average path gains of cluster 2, in dB
    pdb2 = [-inf -inf -3.2 -6.3 -9.4 -12.5 -15.6 -18.7 -21.8];
    % Total average path gains for both clusters, in dB
    pdb = 10*log10(10.^(pdb1/10)+10.^(pdb2/10));
    
    fd = 3;             % Maximum Doppler shift for all paths (identical)
    ds = doppler.bell;  % Bell doppler spectrum, with default parameters
    
    
    
    % Element spacing at the transmit and receive antennas (normalized by the
    % wavelength)
    TxSpacing = 0.5;
    RxSpacing = 0.5;
    
    % Spatial parameters on transmitter side:
    %   Angular spreads - Cluster 1
    AS_Tx_C1 = [14.4 14.4 14.4 14.4 14.4 -inf -inf -inf -inf];
    %   Angular spreads - Cluster 2
    AS_Tx_C2 = [-inf -inf 25.4 25.4 25.4 25.4 25.4 25.4 25.4];
    %   Mean angles of departure - Cluster 1
    AoD_C1 = [225.1 225.1 225.1 225.1 225.1 -inf -inf -inf -inf];
    %   Mean angles of departure - Cluster 2
    AoD_C2 = [-inf -inf 106.5 106.5 106.5 106.5 106.5 106.5 106.5];
    
    % Spatial parameters on receiver side:
    %   Angular spreads - Cluster 1
    AS_Rx_C1 = [14.4 14.4 14.4 14.4 14.4 -inf -inf -inf -inf];
    %   Angular spreads - Cluster 2
    AS_Rx_C2 = [-inf -inf 25.2 25.2 25.2 25.2 25.2 25.2 25.2];
    %   Mean angles of arrival - Cluster 1
    AoA_C1 = [4.3 4.3 4.3 4.3 4.3 -inf -inf -inf -inf];
    %   Mean angles of arrival - Cluster 2
    AoA_C2 = [-inf -inf 118.4 118.4 118.4 118.4 118.4 118.4 118.4];
    
    % Calculation of transmit and receive correlation arrays
    [TxCorrelationMatrix, RxCorrelationMatrix] = ...
        calculateCorrMatrix(Nt, Nr, pdb1, pdb2, TxSpacing, RxSpacing, ...
        AS_Tx_C1, AS_Tx_C2, AoD_C1, AoD_C2, ...
        AS_Rx_C1, AS_Rx_C2, AoA_C1, AoA_C2);
    
    chan = comm.MIMOChannel( ...
        'SampleRate',                Rs, ...
        'PathDelays',                tau, ...
        'AveragePathGains',          pdb, ...
        'MaximumDopplerShift',       fd, ...
        'DopplerSpectrum',           ds, ...
        'TransmitCorrelationMatrix', TxCorrelationMatrix, ...
        'ReceiveCorrelationMatrix',  RxCorrelationMatrix, ...
        'RandomStream',              'mt19937ar with seed', ...
        'Seed',                      99, ...
        'PathGainsOutputPort',       false);
    
     y = step(chan,input.');
    output=awgn(y.',SNR,'measured');
    
elseif strcmp('flat_fading_mimo',channel)   
    % Rs: Input Sample Rate
    Rs = 20e6;
    % Nt: Number of transmit antennas
    Nt=2;
    % Nr: Number of receive antennas
    Nr=2;
    
    chan = comm.MIMOChannel('SampleRate', Rs, ...
            'MaximumDopplerShift', 0, ...
            'PathDelays', 0,...
            'AveragePathGains', 0,...
            'RandomStream', 'mt19937ar with seed',...
            'Seed', 99,...
            'NumTransmitAntennas', Nt,...
            'TransmitCorrelationMatrix', eye(Nt),...
            'NumReceiveAntennas', Nr,...
            'ReceiveCorrelationMatrix', eye(Nr),...
            'PathGainsOutputPort', false,...
            'NormalizePathGains', true,...
            'NormalizeChannelOutputs', true);
    
    y = step(chan,input.');
    output=awgn(y.',SNR,'measured');
else
    output = input;
end
    
    
end

