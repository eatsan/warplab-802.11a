%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Laboratory Algorithmic Research in Network Information (ARNI), 2015
% AUTHORS: Emre ATSAN
%
% See the LICENSE.TXT file for more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef PHY_OFDM_MIMO < IEEE_802_11_PHY
    %PHY_OFDM_MIMO IEEE 802.11 OFDM PHY implementation class with MIMO
    %extension
    %   Inherits IEEE_802_11_PHY abstract class.
    %   This class should contain a concrete class implementation for PLCP,
    %   PMD sublayers (i.e. PLCP_OFDM_MIMO, PMD_OFDM_MIMO).
    
    properties
        SysVars % General system variables for the PHY layer.
        PLCP    % PLCP layer implementation
        PMD     % PMD layer implementation
        PHY_params % PHY layer parameters
    end
    
    methods
        function phy_ofdm_m = PHY_OFDM_MIMO(sys_var, phy_p, tx_id, Nt, Nr)
            
            %Create a new OFDM-MIMO PLCP layer 
            phy_ofdm_m.PLCP = PLCP_OFDM_MIMO(sys_var, tx_id, Nt, Nr);
            
            % Depending on the hardware used, choose the correct PMD layer
            % implementation to instantiate. 
            switch sys_var.rf_device
                case 'sim'
                    phy_ofdm_m.PMD = PMD_OFDM_sim(sys_var);
                case 'warpv2'
                    phy_ofdm_m.PMD = PMD_OFDM_WARPv2(sys_var);
                case 'warpv3'
                    %phy_ofdm.PMD = PMD_OFDM_warpv3(sys_var);
                    error('WARPv3 RF is to be implemented!');
                otherwise
                    error('Unexpected RF_DEVICE type!!');
            end
            
            phy_ofdm_m.PHY_params = phy_p;
            
            %Set general system variables for this PHY object (i.e. CH_BW)
            phy_ofdm_m.SysVars = sys_var;
            
            
        end
        
    end
    
end


