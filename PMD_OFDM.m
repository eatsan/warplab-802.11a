%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Laboratory Algorithmic Research in Network Information (ARNI), 2015
% AUTHORS: Emre ATSAN
%
% See the LICENSE.TXT file for more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef PMD_OFDM < PMD
   %PMD_OFDM OFDM PHY PMD implementation class
   % 
   %   Inherits PMD abstract class. 
    
    properties
        SysVar
    end
    
    methods
        
        %Constructor
        function pmd=PMD_OFDM(sysvar)
            
            pmd.SysVar = sysvar;
            
        end

        
    end
    
end

