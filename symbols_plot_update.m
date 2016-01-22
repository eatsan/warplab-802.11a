%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Laboratory Algorithmic Research in Network Information (ARNI), 2015
% AUTHORS: Emre ATSAN
%
% See the LICENSE.TXT file for more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function symbols_plot_update(rx_syms, hdeMod)
%symbols_plot_update Function which plots & updates input QAM symbols
%
global hScope;

if isempty(hScope)
    %  Scale the constellation for average power equal to 1.
    const_scale  = modnorm(hdeMod.Constellation,'avpow',1); 
    hScope= commscope.ScatterPlot('SamplesPerSymbol',1,...
        'Constellation',const_scale*hdeMod.Constellation, 'RefreshPlot','on');
    hScope.PlotSettings.Constellation = 'on';
    hScope.PlotSettings.ConstellationStyle = 'rd';
    
    update(hScope, rx_syms);
    autoscale(hScope);
    
else
        reset(hScope);
        update(hScope, rx_syms);
        autoscale(hScope);
end

end
