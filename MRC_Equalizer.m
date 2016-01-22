%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Laboratory Algorithmic Research in Network Information (ARNI), 2015
% AUTHORS: Emre ATSAN
%
% See the LICENSE.TXT file for more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ x_estimate ] = MRC_Equalizer( y, h)
%MRC_Equalizer Maximal Ratio Combining equalizer (vector) for given received vector
%and channel estimation vector. Receive diversity gain is maximal.
% Inputs:
% - y = received symbols vector per antenna. Row vector.
% - h = channel estimate vector per antenna. Column vector.


if(size(h,1) > 1) % it is not a column vector, reshape it.
    h = reshape(h,1,size(h,1));
end

x_estimate = (conj(h)*y)./sum(abs(h).^2,2);

end

