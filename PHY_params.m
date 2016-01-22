%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Laboratory Algorithmic Research in Network Information (ARNI), 2015
% AUTHORS: Emre ATSAN
%
% See the LICENSE.TXT file for more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef PHY_params
    %PHY_PARAMS PHY layer specific parameters superclass 
    %   This class is inherited for concrete PHY layer specific parameters
    %   class file is generated. It also provided some common attributes
    %   for the different PHY_params class implementations.
    
    properties
        PreambleDuration; %in seconds
        PLCPHeaderDuration; %in seconds
        MaxMACFrameLength; %in bytes(octats).
    end
    
    methods
        %Class constructor
        function phy_p=PHY_params(preamble_dur, PLCPhdr_dur, maxframe_length)
            if (nargin == 3)
                phy_p.PreambleDuration=preamble_dur;
                phy_p.PLCPHeaderDuration=PLCPhdr_dur;
                phy_p.MaxMACFrameLength=maxframe_length;
            else
                err('PHY_params: Constructor does not have enough input args.');
            end
        end
    end
    
end %PHY_params

