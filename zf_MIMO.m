%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Laboratory Algorithmic Research in Network Information (ARNI), 2015
% AUTHORS: Emre ATSAN
%
% See the LICENSE.TXT file for more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [ IQNorm_RX_MIMO ] = zf_MIMO( nTXNodes, nsym_OFDM, nsubc_data, CE_TX, Data, h_demod, constellation_scale )
%ZF_MIMO Zero-forcing MIMO decoder with successive interference cancellation.
% Ordering of symbol cancellation is done per stream SNR basis. 

IQNorm_RX_MIMO = zeros(nsubc_data,nsym_OFDM,nTXNodes);
%SNR_Estimate_perdatasubc_perOFDMSymbol_perTXdata = zeros(nsubc_data,nsym_OFDM,nTXNodes);
W = zeros(nTXNodes,nTXNodes,nsym_OFDM*nsubc_data);
H = zeros(nTXNodes,nTXNodes,nsym_OFDM*nsubc_data);
y = zeros(nTXNodes,1,nsym_OFDM*nsubc_data);
%snr_sum_rx_antennas = zeros(nsubc_data);
h_mod =  modem.qammod(h_demod);

%CM = zeros(nsubc_data,nsym_OFDM);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MIMO ZF-SIC decoder. (with MRC combining after SIC)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ofdm_symbol=1:nsym_OFDM
        for subc=1:nsubc_data
            linear_seq_id = ((ofdm_symbol-1)*nsubc_data)+subc;

            for stream=1:nTXNodes
                for antenna=1:nTXNodes
                    H(antenna,stream,linear_seq_id) = CE_TX(subc,stream,antenna);
                end
                
                y(stream,1,linear_seq_id) = Data(subc,ofdm_symbol,stream); 
            end

            %Condition number check to validate enough diversity / linear
            %independence exists for the channel matrix observed.
%             CM(subc,ofdm_symbol) = cond(H(1:nTXNodes,1:nTXNodes,linear_seq_id)); 
%             ch_n = mag2db(abs(CM(subc,ofdm_symbol)));
%             
%             if(ch_n > 15)
%                 fprintf('SUBC:[%d]/SYM:[%d] -cond(H): %3f dB. \n',subc,ofdm_symbol,ch_n);
%             end
            
            % ZF equalizer matrix, W
            if (nTXNodes==2)
                W(1:nTXNodes,1:nTXNodes,linear_seq_id) =  fast_2by2_inverse(H(1:nTXNodes,1:nTXNodes,linear_seq_id)); %pinv(H); %pseudo-inverse of channel matrix for subcarrier,subc.
            else
                W(1:nTXNodes,1:nTXNodes,linear_seq_id) = pinv(H(1:nTXNodes,1:nTXNodes,linear_seq_id));
            end
            
%             if(ofdm_symbol==1)
%                 snr_sum_rx_antennas(subc) = sum(SNR(subc,1:nTXNodes));
%             end
        end
        
    end
    
    % Succesive Inteference Cancellation with Optimal Ordering 
    %
    % 1) based on post detection SNR [from book MIMO-OFDM Wireless Communications
    % with MATLAB, ebooks folder]. (SNR is computed by computing the
    % 1/Norm(rows of W).
    %
    % 2) Column norm ordering detection
    for stream=1:nTXNodes
        %1) 
        % norm_W(stream,:) = sqrt( sum(abs(squeeze(W(stream,1:nTXNodes,:))).^2,1));
        
        %2)
         norm_W(stream,:) = sqrt(sum(abs(squeeze(H(:,stream,:))).^2,1));
        
    end
    
    [~,estimation_order]=sort(norm_W,1,'descend');
    %estimation_order(end:-1:1,:);

    x_hat=zeros(nTXNodes,nsym_OFDM*nsubc_data);
    %snr_zf=zeros(nTXNodes,nsym_OFDM*nsubc_data);
  
    for s = 1:nTXNodes
        
         % Spatial stream is detected due to its post-detection SNR (1/W(row)) order or column norm (norm(H(:,col)))detection order.  
        next_stream(:,1) = squeeze(estimation_order(s,:)).';
        
        %If this is the last pass to estimate the next stream per subcarrier, use MRC.
        if(s==nTXNodes)
            
            for i=1:1:nsym_OFDM*nsubc_data
                H_tilde = H(:,next_stream(i),i);
                % Use MRC to combine the estimated version of the last
                % stream from all antennas.
                x_hat(next_stream(i),i) = MRC_Equalizer(y(1:nTXNodes,1,i),H_tilde);
                
                % [sbc,~]=ind2sub([nsubc_data,nsym_OFDM],i); % SLOW, NOT
                % BUILT-IN FUNCTION.
                %sbc = rem(i-1,nsubc_data)+1;
                %snr_zf(next_stream(i),i) = (snr_sum_rx_antennas(sbc).*(1/nTXNodes)) * (1/((norm_W(next_stream(i),i)).^2));
            end
            
            idx=sub2ind(size(x_hat),next_stream',1:1:(nsubc_data*nsym_OFDM));
            x_sliced = QAM_Slicer(h_mod,h_demod,constellation_scale,x_hat(idx));
            
            if(isreal(x_sliced)) %In case of BPSK, output of slicer is not complex number.
                x_sliced = complex(x_sliced);
            end
            
            %Check the EVM value of the first estimated streams per subcarrier.
            rms_evm=step(comm.EVM,squeeze(x_sliced).',squeeze(x_hat(idx)).');
            fprintf('EVM of 2nd estimated stream per subcarriers: %4f percent.\n', rms_evm);
            
        else
            
            for i=1:1:nsym_OFDM*nsubc_data
               x_hat(next_stream(i),i)= W(next_stream(i),1:nTXNodes,i) * y(1:nTXNodes,1,i);
                %[sbc,~]=ind2sub([nsubc_data,nsym_OFDM],i); % SLOW- NOT
                %BUILT-IN
                %sbc = rem(i-1,nsubc_data)+1;
                %snr_zf(next_stream(i),i) = (snr_sum_rx_antennas(sbc).*(1/nTXNodes)) * (1/((norm_W(next_stream(i),i)).^2));
                
            end
            
            idx=sub2ind(size(x_hat),next_stream',1:1:(nsubc_data*nsym_OFDM));
            x_sliced = QAM_Slicer(h_mod,h_demod,constellation_scale,x_hat(idx));
            
            %Check the EVM value of the first estimated subcarriers set.
            
            if(isreal(x_sliced)) %In case of BPSK, output of slicer is not complex number.
                x_sliced = complex(x_sliced);
            end
            
            rms_evm=step(comm.EVM,squeeze(x_sliced).',squeeze(x_hat(idx)).');
            fprintf('EVM of 1st estimated stream per subcarriers: %4f percent.\n', rms_evm);
            
            %Make an hard-decision for the symbols of the first estimation
            %if rms_evm is below 15% %This should be verified.
%             if(rms_evm < 15)
%                 x_hat(idx) = x_sliced;
%             end
            
             for i=1:1:nsym_OFDM*nsubc_data
                y(:,:,i) = y(:,:,i)-H(:,next_stream(i),i)*x_sliced(1,i); % Interference subtraction - SIC
             end

        end
    end
    
    for i=1:1:nsym_OFDM*nsubc_data
        % [sbc,sym]=ind2sub([nsubc_data,nsym_OFDM],i); %SLOW - NOT BUILT IN
        sbc = rem(i-1,nsubc_data)+1;
        sym = (i-sbc)/nsubc_data + 1;
        IQNorm_RX_MIMO(sbc,sym,1:nTXNodes) = x_hat(:,i);
        %SNR_Estimate_perdatasubc_perOFDMSymbol_perTXdata(sbc,sym,1:nTXNodes)=snr_zf(:,i);
        
    end






end

