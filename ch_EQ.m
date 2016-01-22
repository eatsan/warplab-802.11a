%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Laboratory Algorithmic Research in Network Information (ARNI), 2015
% AUTHORS: Emre ATSAN
%
% See the LICENSE.TXT file for more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function signal_equalized = ch_EQ (signal, h_est)
    % Channel Equalizer used for channel compensation.
    % Input: signal = [NSubcarrier x NOFDM_symbols x N_rx_chains] sized matrix after FFT
    %        h_est  = NOFDM_symbols x sized vector for subcarrier based channel est. 
    
    if (size(signal, 1) ~= size(h_est,1))
        error('ch_EQ: signal and channel estimate dimension mismatch!');
    end
    nData_subc = size(signal,1);
    nOFDM_sym = size(signal, 2);
    
    nRX_chains = size(signal, 3);
    
    y= reshape(signal,nData_subc*nOFDM_sym,nRX_chains);
    h = squeeze(h_est);
    signal_equalized=MRC_Equalizer( y, h.');


end