%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Laboratory Algorithmic Research in Network Information (ARNI), 2015
% AUTHORS: Emre ATSAN
%
% See the LICENSE.TXT file for more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef PLCP_ErrorStats < handle
    %ERRORSTATS Error statistics object 
    %   Can be handed in from layer to layer to keep a statistics of errors
    %   over the course of simulations or tests. 
    
    properties
        Detection_error %Frame not detected from LTS and STS.
        STS_Detection_error %Frame not detected from Short Training Sequence
        SIGNAL_not_valid_error % Frame detected, SIGNAL field decoding failed. 
    end
    
    methods
         % Constructor
         function es = PLCP_ErrorStats()
            es.Detection_error=0;
            es.STS_Detection_error =0;
            es.SIGNAL_not_valid_error=0;
         end
         
         function add_DetectionErr(this)
            this.Detection_error=this.Detection_error+1;
         end
         
         function add_STSDetectionErr(this)
             this.STS_Detection_error=this.STS_Detection_error+1;
         end
         
         function add_SIGNALNVErr(this)
             this.SIGNAL_not_valid_error=this.SIGNAL_not_valid_error+1;
         end

    end
    
end

