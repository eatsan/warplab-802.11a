%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Laboratory Algorithmic Research in Network Information (ARNI), 2015
% AUTHORS: Emre ATSAN
%
% See the LICENSE.TXT file for more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef (Abstract) PLCP < handle
    %PLCP Physical Layer Convergence Procedure(PLCP) sublayer abstract class
    %   Every PHY should provide an implementation for this sublayer abstract class
    
    properties
    end
    
    methods (Abstract, Static)
        
        PLCP_TX
        PLCP_RX

    end
    
end

