%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Laboratory Algorithmic Research in Network Information (ARNI), 2015
% AUTHORS: Emre ATSAN
%
% See the LICENSE.TXT file for more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef (Abstract) TXVECTOR < handle
    %TXVECTOR TXVECTOR parameters for PLCP & MAC layer management
    %
    %    TXVECTOR supplies the PHY with per-packet transmit parameters.
    %    This is an abstract class to be implemented depending on the PHY
    %    implementation (OFDM-PHY, HT-PHY, etc.)
    
    properties
    end
    
    methods
    end
    
end

