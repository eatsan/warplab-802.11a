%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Laboratory Algorithmic Research in Network Information (ARNI), 2015
% AUTHORS: Emre ATSAN
%
% See the LICENSE.TXT file for more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ invM ] = fast_2by2_inverse( M )
%FAST_2by2_INVERSE For the special case of 2 by 2 matrices, finding the
%inverse can be implemented faster in the following way: 
invM =  ( 1 / (M(1,1)*M(2,2) - M(1,2)*M(2,1))) .* [M(2,2), -M(1,2) ; -M(2,1), M(1,1)  ];


end

