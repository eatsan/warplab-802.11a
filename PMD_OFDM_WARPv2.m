%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Laboratory Algorithmic Research in Network Information (ARNI), 2015
% AUTHORS: Emre ATSAN
%
% See the LICENSE.TXT file for more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef PMD_OFDM_WARPv2 < PMD_OFDM
    %PMD_OFDM_WARPv2 Warp v2 radio PMD layer
    % Interpolates the 20Mhz signal to 40Mhz, reception, transmission logic
    % with the hardware
    %  
    
    properties
        INTERP_RATE = 2;        %Interpolation rate (1 or 2)
        TX_SCALE = 1.0;         %Scale for Tx waveform ([0:1])
        interp_filt2;           %Interpolation filter
        TX_NUM_SAMPS;           %Number of samples for I/Q queue of WARPv2 hardware.
    end
    
    methods
         %Constructor
        function pmd=PMD_OFDM_WARPv2(sysvar)
            
            %Call superclass constructor
            pmd = pmd@PMD_OFDM(sysvar);
            
            pmd.TX_NUM_SAMPS = sysvar.IQ_queue_size_samples;
            %Define a halfband 2x interp filter response
            pmd.interp_filt2 = zeros(1,43);
            pmd.interp_filt2([1 3 5 7 9 11 13 15 17 19 21]) = [12 -32 72 -140 252 -422 682 -1086 1778 -3284 10364];
            pmd.interp_filt2([23 25 27 29 31 33 35 37 39 41 43]) = pmd.interp_filt2(fliplr([1 3 5 7 9 11 13 15 17 19 21]));
            pmd.interp_filt2(22) = 16384;
            pmd.interp_filt2 = pmd.interp_filt2./max(abs(pmd.interp_filt2));
            
        end
        
        
        % PLCP RX procedure
        function [PPDU] = PMD_RX(this, rx_IQ)
            dim = size(rx_IQ,2);
            for i=1:dim
                rx_IQ(:,i)=warplab_correctdco(rx_IQ(:,i),480);
            end
            DECIMATE_RATE = this.INTERP_RATE;
            % Decimate
            if(DECIMATE_RATE == 1)
                raw_rx_dec = rx_IQ;
            elseif(DECIMATE_RATE == 2)
                raw_rx = filter(this.interp_filt2, 1, rx_IQ);
                raw_rx_dec(:,:) = raw_rx(1:2:end,:);
            end
            
            PPDU = raw_rx_dec; 
        end
        
        % PLCP TX procedure
        function tx_vec_samples = PMD_TX(this, PPDU)
            
            %Pad with zeros for transmission
            tx_vec_padded = [PPDU zeros(1,(this.TX_NUM_SAMPS/this.INTERP_RATE)-length(PPDU))];
            
            % Interpolate
            if(this.INTERP_RATE == 1)
                tx_vec_samples = tx_vec_padded;
            elseif(this.INTERP_RATE == 2)
                tx_vec_2x = zeros(1, 2*numel(tx_vec_padded));
                tx_vec_2x(1:2:end) = tx_vec_padded;
                tx_vec_samples = filter(this.interp_filt2, 1, tx_vec_2x);
            end
            
            %Scale the Tx vector
            tx_vec_samples = this.TX_SCALE .* tx_vec_samples ./ max(abs(tx_vec_samples));
            
        end
        
    end
    
end

