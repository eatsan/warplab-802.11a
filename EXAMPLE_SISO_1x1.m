%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% � All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Laboratory Algorithmic Research in Network Information (ARNI), 2015
% AUTHORS: Emre ATSAN
%
% See the LICENSE.TXT file for more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;

%-------------------------------------------------------------------------
%
% This example uses WARPLab-wifi package for generating digital BB samples
% and WARPLab Reference Design (>v7.4.0) for transferring samples to WARPv2
% radios. 
%
% Please see WARP Project website for how to setup WARPLab Framework. 
% http://warpproject.org/trac/wiki/WARPLab/QuickStart
%
%-------------------------------------------------------------------------

% ------------------------------------------------------------------------
% Declare system setup variables
% ------------------------------------------------------------------------
% Number of transmitting antennas 
sysvar.N_Tx=1; % This is a Single-Input Single-Output example 
% Number of receiving antennas
sysvar.N_Rx=1;
% Channel Bandwidth (only 20 MHz supported at this point)
sysvar.Ch_BW=20; %MHz
% Guard Interval (GI) option for OFDM operation (Yet no support for short GI)
sysvar.short_GI=0; % 0: Long GI, 1: Short GI
% Maximum MSDU length (i.e., maximum payload length) 
sysvar.maxLength = 4095; %bytes
% Decision Statistics for detection threshold from LTS training
sysvar.LTSCorr_THRESH = 0.8; % Keep this value, if you are not sure what it is. 
% Number of Cyclic Prefix (CP) samples to use in FFT (for OFDM operation)
sysvar.FFT_OFFSET = 4; % Keep this value, if you are not sure what it is.
% Collect PCLP Error statistics during the test
sysvar.PLCPErrorStats = true; 
% Maximum number of samples that we can read/write from the radio I/Q stack.
sysvar.IQ_queue_size_samples = 2^14; % For WARPv2 2^14 samples.
% Hardware for RF transmission (Yet no support for warpv3.)
sysvar.rf_device = 'warpv2'; 
% Plot DATA field symbols before channel decoding.
% Use this variable to see received QAM symbols on the transmit
% constellation. Very useful for debugging.
sysvar.plot_data_symbols=true;

% ------------------------------------------------------------------------
% Declare WARPLab related variables
% ------------------------------------------------------------------------
% Number of WARP nodes used in this example
NUMNODES = 2; % 1 node for source, 1 node for destination. 
% Use Automatic Gain Control (AGC) for RX gains?
USE_AGC = true;
% Operating Frequency Band 
FREQ_BAND = 2.4;
% ISM Band channel to use
CHANNEL = 4; 

% ------------------------------------------------------------------------
% Cyclic Redundancy Check (CRC) variables
% ------------------------------------------------------------------------
% In this example, we will check the integrity of the transmitted
% information using CRC-32 error detection. This is a standard operation in
% 802.11a OFDM PHY (2012). 

% ITU-T CRC 32 for FCS settings
% Generator polynomial defined in the Std. 802.11-2012 - page 400
CRC32_PolyCoeff = zeros(1,33);
CRC32_PolyCoeff([32 26 23 22 16 12 11 10 8 7 5 4 2 1 0]+1) = 1;
CRC32_PolyCoeff = fliplr(CRC32_PolyCoeff);

% Create a CRC generator and a CRC detector object.
h_crcgenerator = crc.generator('Polynomial', CRC32_PolyCoeff, ...
    'InitialState', ones(1,32), ...
    'FinalXOR', ones(1,32));
h_crcdetector = crc.detector('Polynomial', CRC32_PolyCoeff, ...
    'InitialState', ones(1,32), ...
    'FinalXOR', ones(1,32));

% ------------------------------------------------------------------------
%  WARPLab Framework initialization
% ------------------------------------------------------------------------
%  See WARPLab documentation and examples for the details of these
%  operation. Basically WARPLab lets us to connect and transmit/receive raw
%  samples to WARP node radio queues using MATLAB commands.

%Create a vector of node objects
nodes = wl_initNodes(NUMNODES);

%Create a UDP broadcast trigger and tell each node to be ready for it
eth_trig = wl_trigger_eth_udp_broadcast;
wl_triggerManagerCmd(nodes,'add_ethernet_trigger',[eth_trig]);

% Get IDs for the radio interfaces on the boards. 
% In this example we assume that each WARP node has 2 radio cards
% installed. 
[RFA1, RFB1] = wl_getInterfaceIDs(nodes(1)); % Node 1
[RFA2, RFB2] = wl_getInterfaceIDs(nodes(2)); % Node 2

% Set Node 1: transmitter, and Node 2: receiver. 
% Make sure that WARP radios have correct jumper settings for Node 1 and
% Node 2. (See WARP documentation for that.)
if (NUMNODES >= 2)
    node_tx = nodes(1);
    node_rx = nodes(2);
else
    error('ERROR: Need at least two WARP nodes!');
end
    

% Use RFA1 for Node 1 tranmission. 
RF_TX = RFA1; % 1x1 SISO 
% Use RFA2 for Node 2 reception.
RF_RX = RFA2; 
% Note: Make sure that the RF card which is NOT used is ALSO connected to
% a 50 Ohm terminator or antenna.

% Set the channel and frequency band information for all RF interfaces.
wl_interfaceCmd(nodes,'RF_ALL','channel',FREQ_BAND,CHANNEL);

% Set the transmit power (See WARPLab documentation.)
TX_BB_GAIN = 0;  % Baseband gain setting
TX_RF_GAIN = 15; % RF gain setting

%Compute estimated RF output power in dBm. (MAX2829 documentation).
if(FREQ_BAND==2.4)
    est_tx_power_dBm = -6 + (TX_BB_GAIN*1.5-4.5)+(TX_RF_GAIN*0.5-32)+30;
elseif(FREQ_BAND==5)
    est_tx_power_dBm = -6 + (TX_BB_GAIN*1.5-4.5)+(TX_RF_GAIN*0.5-34)+30;
end
wl_interfaceCmd(node_tx,'RF_ALL','tx_gains',TX_BB_GAIN,TX_RF_GAIN);

% Set the Automatic Gain Control (AGC) for radios
if(USE_AGC)
    % Refer to WARP Project documentation for the details.
    
    % We are using AGC
    wl_interfaceCmd(node_rx,'RF_ALL','rx_gain_mode','automatic');
   
    wl_basebandCmd(node_rx,'agc_target',-11);
    % Let's give a bit of delay before triggering AGC. 
    % This is because TX & RX nodes might not be completely synchronized.
    wl_basebandCmd(node_rx,'agc_trig_delay', 150);
    % DCO Correction
    wl_basebandCmd(node_rx,'agc_dco', true);
    % RX/TX low pass filter settings at the RF. (Refer to documentation for
    % the values. No need to change at this point). 
    wl_interfaceCmd(nodes,'RF_ALL','rx_lpf_corn_freq',2);
    wl_interfaceCmd(nodes,'RF_ALL','tx_lpf_corn_freq',1);
    % Use this to add delay manually to the transmitter after Ethernet
    % trigger.
    wl_basebandCmd(nodes,'tx_delay',0); % No delay
    
    % Read maximum # of TX samples from the WARP radios. 
    sysvar.IQ_queue_size_samples = node_tx.baseband.txIQLen;
    wl_basebandCmd(nodes,'tx_length',sysvar.IQ_queue_size_samples);

    % Sampling Frequency (For WARPv2 40 MHz).
    sysvar.SAMP_FREQ = wl_basebandCmd(node_tx,'tx_buff_clk_freq');
else
    % You can manually set the RX Gain values if necessary. Refer to
    % WARPLab documentation for it. We will use AGC in this example.
    error('Please USE_AGC=true for this example.');
end


% ------------------------------------------------------------------------
%  TX Node SETUP
% ------------------------------------------------------------------------

% Instantiate a PHY_params object to store Physical layer variables
% We only need maxLength at this point. It can be extended later.
phy_p = PHY_params(0, 0, sysvar.maxLength); 


% Create a OFDM PHY Layer object to use for the transmitter.
phy_tx = PHY_OFDM(sysvar,phy_p);

% Length of the MAC layer payload to be generated. 
tx_length = 250;  %bytes
% 802.11a PHY transmission rate (Choose from 6,9,12,18,24,36,48,54).
tx_datarate = 24; %Mbps
% TX Vector (from 802.11 OFDM PHY spec.) 
txv = TXVECTOR_OFDMPHY(tx_length, tx_datarate, phy_p);

% ------------------------------------------------------------------------
%  RX Node SETUP
% ------------------------------------------------------------------------

% Creating a OFDM PHY Layer object to use for the receiver.
phy_receiver = PHY_OFDM(sysvar,phy_p);


% ------------------------------------------------------------------------
%  OVER-THE-AIR Operation
% ------------------------------------------------------------------------

%  --- TX NODE Operation

% Initialize MAC Protocol Data Unit (MPDU) vector. 
MPDU_T = zeros(txv.LENGTH,1); % our payload

% Fill with random bits for test purposes.
MPDU_T(1:txv.LENGTH-32,1) = randi([0 1],txv.LENGTH-32,1);

% Concatenate CRC-32 at the end of MPDU
MPDU_T(1:txv.LENGTH,1) = generate(h_crcgenerator, MPDU_T(1:txv.LENGTH-32,1));
 
% Deliver TXVECTOR and payload to TX procedure of PLCP layer.
PPDU = phy_tx.PLCP.PLCP_TX(txv, MPDU_T); 

% Create baseband signal for RF transmission (WARP)
TX_IQ(1,:) = phy_tx.PMD.PMD_TX(PPDU);


%Write the Tx waveform to the Tx warp node
wl_basebandCmd(node_tx,RF_TX, 'write_IQ', TX_IQ.');

%Enable the Tx and Rx radios
wl_interfaceCmd(node_tx,RF_TX,'tx_en');
wl_interfaceCmd(node_rx,RF_RX,'rx_en');

%Enable the Tx and Rx buffers
wl_basebandCmd(node_tx,RF_TX,'tx_buff_en');
wl_basebandCmd(node_rx,RF_RX,'rx_buff_en');

%Trigger the Tx/Rx cycle at both nodes
eth_trig.send();
    
%  --- RX NODE Operation

% Retrieve the received waveform from the Rx node
RX_IQ = wl_basebandCmd(node_rx,RF_RX,'read_IQ', 0, sysvar.IQ_queue_size_samples);
rx_RSSI = wl_basebandCmd(node_rx,RF_RX,'read_RSSI');

% Disable the Tx/Rx radios and buffers
wl_basebandCmd(nodes,'RF_ALL','tx_rx_buff_dis');
wl_interfaceCmd(nodes,'RF_ALL','tx_rx_dis');

% Received signal decoding
[PPDU_R] = phy_receiver.PMD.PMD_RX(RX_IQ);
[rxv, MPDU_R] = phy_receiver.PLCP.PLCP_RX(PPDU_R);

% ------------------------------------------------------------------------
%  EVALUATION 
% ------------------------------------------------------------------------

% Check if the received bits are correct?
err = false; 

if(any(MPDU_R ~=-1)) % Any error in signal detection or wrong SIGNAL field 
    % in PLCP_RX() returns -1 value for MPDU_R vector.
    %Check CRC-32
    crc_err = 0; % Initialize to false.
    [octets, crc_err] = detect(h_crcdetector, double(MPDU_R));
    
    if(crc_err == 1)
        err=true;
    end
 
else
    err = true; 
end

% Print the estimated TX power value. (depends on TX, RF gains).
fprintf('[TX_POWER] Average transmit power: %d dBm.\n',est_tx_power_dBm);

% Print MPDU (Payload) information
fprintf('[DATARATE] RATE: %d Mbps.\n',tx_datarate);
fprintf('[FRAME_LENGTH] SIZE: %d bytes.\n',tx_length);

% Print the result.
if (err)
    fprintf('[RESULT] Received frame dropped!\n');
else
    fprintf('[RESULT] Received frame successfully decoded!\n');
end

% OUTPUT the PLCP Error Statistics.
% phy_receiver.PLCP.ErrStats

% Plot PSD of the received signal.
figure;
% Ignore the first 400 samples since it might contained before AGC samples
[Pxx,W] = pwelch(RX_IQ(400:length(PPDU)*2,1),[],[],4096,40);
%calculates the power density function(Pxx) of the signal to be transmitted and the frequency vector W
plot([-2048:2047]*40/4096,10*log10(fftshift(Pxx)),'b');
%plots the PSD of signal in dB after shifting the zero frequency component to the center.
xlabel('Baseband Frequency (MHz)')
ylabel('Power Spectral Density (dB)')
title('Received signal');
grid on

% Compute PAPR (Peak-to-average-power ratio) for the transmitted signal.
abs_tx_iq = abs(TX_IQ(1,1:(length(PPDU)*2)));
fprintf('[TX_PAPR] PAPR = %4f dB \n',10*log10(max(sqrt(abs_tx_iq))/mean(sqrt(abs_tx_iq))));

% Total transmit time (in musec.)
tx_time = (length(PPDU)*5.0000e-08)/(10^-6);
fprintf('Total Data TX time (SISO): %d musec.\n',tx_time);
