%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Laboratory Algorithmic Research in Network Information (ARNI), 2015
% AUTHORS: Emre ATSAN
%
% See the LICENSE.TXT file for more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef PMD_OFDM_sim < PMD_OFDM
    %PMD_OFDM_SIM Dummy class to test simulation of the PMD layer 
    %  
    
    properties
    end
    
    methods
 
        %Constructor
        function pmd=PMD_OFDM_sim(sysvar)
            
            %Call superclass constructor
            pmd = pmd@PMD_OFDM(sysvar);
            
        end
        
        % PLCP RX procedure
        function [RSSI, PPDU] = PMD_RX(this, rx_IQ)
            %TODO
            PPDU = rx_IQ; 
            RSSI = 1; % Dummy RSSI for simulations.
        end
        
        % PLCP TX procedure
        function tx_samples = PMD_TX(this, PPDU)
            tx_samples = [PPDU(1:end) complex(zeros(1,50),zeros(1,50))];
        end
    end
    
end

