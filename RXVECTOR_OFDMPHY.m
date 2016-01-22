%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Laboratory Algorithmic Research in Network Information (ARNI), 2015
% AUTHORS: Emre ATSAN
%
% See the LICENSE.TXT file for more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef RXVECTOR_OFDMPHY < RXVECTOR
    %RXVECTOR_OFDMPHY RXVECTOR data structure for OFDM-PHY implementation
    %for IEEE Std. 802.11-2012.
    %   This class inherits the abstract RXVECTOR abs. class and implements
    %   the parameter list for the OFDM-PHY (Legacy) case.
    
    properties
        LENGTH
        RSSI
        DATARATE
        SERVICE
        RCPI
        ANT_STATE
        RX_START_OF_FRAME_OFFSET
    end
    
    methods
        % Constructor
        function rv = RXVECTOR_OFDMPHY(length, rate)
            
            %Set PSDU length in bytes.
            rv.LENGTH = length; %bytes
            
            rv.RSSI = 0;%Dummy at this point.
            
            
            %Set Datarate
            if(rate==6 || rate==9 || rate==12 ||rate==18 ||rate==24 || rate==36 ||rate==48 ||rate==54)
                rv.DATARATE = rate;
            else
                error('RXVECTOR_OFDMPHY: DATARATE should be one of the specified in Table 18-1 - Std. IEEE802.11-2012.');
            end
            
            rv.SERVICE = zeros(1,16); % Shall be null
            
            rv.RCPI = 0; %Not available yet.
            
            rv.ANT_STATE = 0; %Not available yet.
            
            rv.RX_START_OF_FRAME_OFFSET=0; %N/A
            %Call superclass constructor (if exists)
        end
    end
    
end

