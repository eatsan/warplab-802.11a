% Copyright (c) 2009, Won Y. Yang All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution
% 
% CHANGES: Edited to support more than one block interleaving
% (Multi-channel support) 
% Copyright (c) 2014, Emre Atsan. All rights reserved.

function x=deinterleaver(y,Nrow,Ncbps,Nbpsc)
% y    : Frame of bits/symbols to be deinterleaved
% Nrow : Number of rows of the block interleaver
% Ncbps: Number of coded bits per symbol
% Nbpsc: Number of coded bits per subcarrier
% x    : Deinterleaved frame of bits/symbols 

y = reshape(y,Ncbps,[]); 

x=y; % To make x have the same size as y
i=0:Ncbps-1; s= max(Nbpsc/2,1); 
j=s*floor(i/s) + mod(i+Ncbps-floor(Nrow*i/Ncbps),s); % 1st permutation
z(j+1,:)=y;  k=i;
i = Ncbps/Nrow*mod(k,Nrow) + floor(k/Nrow); % 2nd permutation
%Ncol=Ncbps/Nrow; i=Nrow*mod(k,Ncol) + floor(k/Ncol); % 2nd permutation
x(i+1,:)=z;

x = reshape(x,[],1); 
end
