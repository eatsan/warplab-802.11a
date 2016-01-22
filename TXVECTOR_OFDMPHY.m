%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Laboratory Algorithmic Research in Network Information (ARNI), 2015
% AUTHORS: Emre ATSAN
%
% See the LICENSE.TXT file for more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef TXVECTOR_OFDMPHY < TXVECTOR
    %TXVECTOR_OFDMPHY TXVECTOR data structure for OFDM-PHY implementation
    %for IEEE Std. 802.11-2012.
    %   This class inherits the abstract TXVECTOR abs. class and implements
    %   the parameter list for the OFDM-PHY (Legacy) case.
    
    properties
        LENGTH
        DATARATE
        SERVICE
        TXPWR_LEVEL
        TIME_OF_DEPARTURE_REQUESTED
    end
    
    methods
        % Constructor
        function tv = TXVECTOR_OFDMPHY(length, rate, phy_params)
            
            %Set PSDU length in bytes.
            if(length <= phy_params.MaxMACFrameLength)
                tv.LENGTH = length * 8; % bytes - > bits
            else
                error('TXVECTOR_OFDMPHY: PSDU size exceeds maximum allowed.');
            end
            
            %Set Datarate
            if(rate==6 || rate==9 || rate==12 ||rate==18 ||rate==24 || rate==36 ||rate==48 ||rate==54)
                tv.DATARATE = rate;
            else
                error('TXVECTOR_OFDMPHY: DATARATE should be one of the specified in Table 18-1 - Std. IEEE802.11-2012.');
            end
            
            tv.SERVICE = zeros(1,16); % 16 bits. All zeros.
            
            tv.TXPWR_LEVEL = 1; %Not available yet.
            
            tv.TIME_OF_DEPARTURE_REQUESTED = false; %Not available yet.
            
            %Call superclass constructor (if exists)
        end
    end
    
end

