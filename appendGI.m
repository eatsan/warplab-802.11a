%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Laboratory Algorithmic Research in Network Information (ARNI), 2015
% AUTHORS: Emre ATSAN
%
% See the LICENSE.TXT file for more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ y ] = appendGI( x, CPLength )
%APPENDGI Adds the cyclic prefix to the beginning of the OFDM symbol
%   Inputs: 
%       * CPLength flag indicates the sample length of the GI.
%       * x: OFDM symbol samples 
%   Outputs:
%       * y: GI appended OFDM symbol 
% 

%Reshape it to column vector
x = reshape(x,[],1);

tx_cp = x((end-CPLength+1 : end),1);

y = [tx_cp;x];

end

