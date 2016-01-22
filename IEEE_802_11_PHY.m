%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Laboratory Algorithmic Research in Network Information (ARNI), 2015
% AUTHORS: Emre ATSAN
%
% See the LICENSE.TXT file for more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef (Abstract) IEEE_802_11_PHY
    %IEEE_802_11_PHY Abstract IEEE 802.11 PHY implementation class
    %
    %   A physical layer (PHY) implementation for IEEE Std 802.11-2012.

    
    properties (Abstract=true)
        PLCP  %-PLCP sublayer abstract class implementation necessary 
              % for a concrete class. (i.e. PLCP_OFDM, PCLP_DSSS) 
        PMD   %-PMD sublayer abstract class implementation necessary 
              % for a concrete class. (i.e. PMD_OFDM, PMD_DSSS)
        PHY_params %PHY parameters superclass for different PHY implementations
    end
    
    methods
    end
    
end

