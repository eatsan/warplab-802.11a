%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Laboratory Algorithmic Research in Network Information (ARNI), 2015
% AUTHORS: Emre ATSAN
%
% See the LICENSE.TXT file for more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [ sliced_symbols ] = QAM_Slicer( h_mod, h_demod, scaleConstellation_avpow1, x_tmp )
% QAM_SLICER Given a complex symbol, this function aims to slice it to the
%  closest QAM symbol from the given modulator.
%  Note: It might not be the fastest slicer in the world. 
%
%  INPUT: 
%        h_mod: modulator object
%        h_demod: demodulator object
%        scaleConstellation_avpow1: scaling factor used for the modulator 
%                                  to make the average power = 1
%        x_temp : received signal that will be sliced. 
%
% OUTPUT:
%        sliced_symbols



% This function basically demodulates the received signal and remodulates
% with the given modulator object
sliced_symbols = modulate(h_mod, ...
                demodulate(h_demod, x_tmp.*(1/scaleConstellation_avpow1))).*scaleConstellation_avpow1; 

end

