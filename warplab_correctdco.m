%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Laboratory Algorithmic Research in Network Information (ARNI), 2015
% AUTHORS: Emre ATSAN, Melissa DUARTE
%
% See the LICENSE.TXT file for more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [RxData_withoutDCO] = warplab_correctdco(RxData_withDCO,AGC_Set_Address)

% DC Offset (DCO) removal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 0: Compute the DC offset in the I and Q paths by averaging these
% signals over the duration of a short symbol. Perform this averaging after
% the AGC has set its gains, the sample at which the AGC is done setting
% the gains is given by AGC_Set_Address.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Average I nand Q data over the duration of a short symbol. Lasts 16
% samples but for WARPLab we are upsampling by 2 so a short symbol lasts 32
% smaples. For WARPLab fs = 40MHz so 32 samples last 0.8us.
RxData_withDCO_I = real(RxData_withDCO);
DCO_I = mean(RxData_withDCO_I(AGC_Set_Address+1:AGC_Set_Address+32));  
RxData_withDCO_Q = imag(RxData_withDCO);
DCO_Q = mean(RxData_withDCO_Q(AGC_Set_Address+1:AGC_Set_Address+32));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Substract the computed I and Q DC offsets from their respective
% signals for the rest of the samples and store result in
% RxData_withoutDCO_1. The first 'AGC_Set_Address'+32 samples remain
% unchanged, these are the samples before the AGC sets the gains and the 32
% samples that are used to compute the DC_Offset in step0 above. We leave
% these samples unchanged so that the user can observe the exact signal 
% output by the ADC's, during the time the AGC is computing the gains, and
% the 32 samples used to estimate the DCO.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute length of vector of received samples
lenRxData = length(RxData_withDCO); 
% The first 'AGC_Set_Address'+32 samples remain unchanged
RxData_withoutDCO_I_1 = RxData_withDCO_I(1:AGC_Set_Address+32);
RxData_withoutDCO_Q_1 = RxData_withDCO_Q(1:AGC_Set_Address+32);
% Remove DCO 
RxData_withoutDCO_I_1(AGC_Set_Address+32+1:lenRxData) = ...
    RxData_withDCO_I(AGC_Set_Address+32+1:lenRxData)-DCO_I;
RxData_withoutDCO_Q_1(AGC_Set_Address+32+1:lenRxData) = ...
    RxData_withDCO_Q(AGC_Set_Address+32+1:lenRxData)-DCO_Q;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Apply a higpass filter to remove any extra DCO. In Step 1 the DCO
% computed in Step 0 was substracted from the received signal. This
% porcessing reduces the DCO but the DCO estimated in Step 0 may not be the
% true DCO, a residual error may remain and it must be removed. This
% residual error is removed applyig a high pass filter.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The highpass filter is designed offline using MATLAB's fdatool. The
% filter coefficients are defined below.
Filter_Num = 0.99843166591671906 * [1,-1];
Filter_Den = [1,-0.99686333183343789];
% Filter the residual DCO
RxData_withoutDCO_I = RxData_withoutDCO_I_1;
RxData_withoutDCO_I(AGC_Set_Address+32+1:lenRxData) = ...
    filter(Filter_Num,Filter_Den,RxData_withoutDCO_I_1(AGC_Set_Address+32+1:lenRxData));
RxData_withoutDCO_Q = RxData_withoutDCO_Q_1;
RxData_withoutDCO_Q(AGC_Set_Address+32+1:lenRxData) = ...
    filter(Filter_Num,Filter_Den,RxData_withoutDCO_Q_1(AGC_Set_Address+32+1:lenRxData));
% Create output vector
RxData_withoutDCO = complex(RxData_withoutDCO_I,RxData_withoutDCO_Q);

