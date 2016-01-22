%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Laboratory Algorithmic Research in Network Information (ARNI), 2015
% AUTHORS: Emre ATSAN
%
% See the LICENSE.TXT file for more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef PLCP_OFDM_MIMO < PLCP_OFDM
    %PLCP_OFDM_MIMO Extended version of PLCP_OFDM for the MIMO transmission
    %
    %   This class aims to extend the legacy OFDM PHY for 802.11a IEEE for
    %   multiple tx chains with minimal changes on the original standard
    %   definitions. Some features of this extension is specified in HT PHY
    %   defined for the 802.11n standard. For the sake of simplicity, we
    %   avoid trying to implement the full technical specification for HT
    %   PHY and use this one for MIMO tx/rx.
    %
    %  NOTE that this class inherits its superclass PLCP_OFDM. This means 
    %  unless it overwrites, it has all the properties and methods of
    %  PLCP_OFDM. 
    %
    properties
        tx_id       % Transmitter ID in the MIMO TX chain.
        cond_n_mimo % Channel Condition number observed in dB
    end
    
    methods (Access = public)
        function p = PLCP_OFDM_MIMO(sysvar, tx_id, Nt, Nr)
            if nargin ~= 4
                error('PLCP_OFDM_MIMO: Not enough input arguments for constructor.');
            end
            
            % Call superclass constructor
            p = p@PLCP_OFDM(sysvar);
            p.tx_id = tx_id;
            
            p.Nt = Nt;
            p.Nr = Nr;
            
            if(p.tx_id > p.SysVar.N_Tx)
                error('PLCP_OFDM_MIMO->PLCP_OFDM_MIMO: Number of TX chains is not set correctly. Check your sysvar.N_Tx variable!');
            end
            
            if(p.SysVar.N_Tx == 2)
                
                % tx_id=0 is assigned to a receiver only PLCP_OFDM_MIMO
                % node.
                if(p.tx_id > 0)
                    %Share the pilot subcarriers
                    p.index_subc_pilots = p.index_subc_pilots(p.tx_id:2:end);
                    p.pilots_f = p.pilots_f(p.tx_id:2:end);
                end

                
                % Give a circular shift to the second transmitter STS
                % field. This is to avoid unexpected signal cancellations
                % in MIMO transmissions. (This is a standard operation in
                % IEEE 802.11 HT OFDM PHY.
                if(p.tx_id ==2)
                    %Circular shift by -200 nsec (-4 samples). Negative
                    %cyclic delay. (CDD - cyclic delay diversity from
                    %802.11 HT PHY.)
                    p.sts_t_ext= circshift(p.sts_t_ext.',-4).';
                    
                end
            
            end
            
            
        end
        
        % PLCP RX procedure
        function [rxvector,PSDU] = PLCP_RX(this, PPDU)
            %Inputs:
            % PPDU: PLCP Protocol Data Unit. Input from PMD.
            %       Have more than 1 vector for # of receive antennas
            %       this.Nr>1.

            
            %Output:
            % PSDU: PLCP Service Data Unit (MAC payload - MPDU).
            % rxvector: RXVECTOR parameters for MAC layer operations.
            
            if(this.Nr ~= size(PPDU,2))
                warning('PLCP_OFDM_MIMO->PLCP_RX: Not enough sample dimensions received from the PMD layer!');
                rxvector=-1; PSDU=-1;
                return;
            end
            
            for rx_antenna=1:1:this.Nr
                % Signal detection from preamble
                % --We add a fudge factor for the search (~100)
                % --because there may be jitters so we search a bit more than
                % --preamble size.
                [decision(rx_antenna), start_of_payload_samples(rx_antenna)] = this.PLCP_Signal_Detection(squeeze(PPDU(1:this.NPreamble_Samples+100,rx_antenna)));
            end
            
            if(all(decision)~=true)
                warning('PLCP_OFDM_MIMO->PLCP_RX: Frame Detection failed!');
                rxvector=-1; PSDU=-1;
                return;
            end
            
            
            for rx_antenna=1:1:this.Nr
                
                %This variable was used to check if it is better to decide
                %using the stronger antenna value for all antennas. However
                %it did not work out, so just ignore it for now. 
                master_antenna =rx_antenna;
                
                %Sanity Check for synchronization
                if ((start_of_payload_samples(master_antenna)-this.NLTS_Samples) < -5) %More than 5 samples off is not a good estimate.
                    warning('PLCP_OFDM_MIMO->PLCP_RX: Frame Detection failed! Wrong start of payload estimation.');
                    rxvector=-1; PSDU=-1;
                    this.ErrStats.add_DetectionErr();
                    return;
                elseif ((start_of_payload_samples(master_antenna)-this.NLTS_Samples) < 1) %Less than 5 samples off can be corrected.
                    start_of_payload_samples(master_antenna) = start_of_payload_samples(master_antenna) + 1 - (start_of_payload_samples(master_antenna)-this.NLTS_Samples);
                end
                
                % Coarse CFO estimation & correction
                PPDU_cfo_cor(:,rx_antenna) = this.PLCP_Coarse_CFO(squeeze(PPDU(:,rx_antenna)),start_of_payload_samples(master_antenna)); %time domain

                %Channel Estimation (using LTS training)
                preamble_lts = PPDU_cfo_cor((start_of_payload_samples(master_antenna)-this.NLTS_Samples):(start_of_payload_samples(master_antenna)-1),rx_antenna);
                h_est(:,:,rx_antenna) = this.PLCP_Channel_Estimate(preamble_lts); % Frequency domain channel estimates
                
                % SNR estimate (using LTS training)
                [SNR_est(:,:,rx_antenna),SNR_est_avg(rx_antenna,:), noise_variance(:,:,rx_antenna) ]= this.PLCP_SNR_est_LTS(preamble_lts);
                
                
            end
            
            %Set the variable to save condition number of the MIMO channel.
            cond_n_dB = zeros(this.NOFDM_data_subc/2,1);
            
            i=0;
            for subc = this.index_subc_data %for all data subcarriers
                i=i+1;
                cond_n_dB(i) = 10*log10(abs(cond(squeeze(h_est(subc, :,:)))));
            end
            
            %Set the condition number of the observed MIMO channel.
            this.cond_n_mimo = cond_n_dB; 
            
            %Set the average SNR per stream per rx antenna.
            this.SNR = SNR_est_avg; 
            
            this.start_of_payload = start_of_payload_samples;
            
            %Extract SIGNAL field & Create RXVECTOR
            [length, d_rate] = this.PLCP_SIGNAL_Decode(PPDU_cfo_cor, h_est, start_of_payload_samples);
            
            if(length==-1 && d_rate==-1)
                warning('PLCP_OFDM->PLCP_RX: SIGNAL not decoded!');
                rxvector=-1; PSDU=-1;
                return;
            else
                rxvector = RXVECTOR_OFDMPHY(length, d_rate);
                % TODO: Ideally this needs to create an event to signal MAC Layer.
                %  -- PHY-RXSTART.indication(RXVECTOR) primitive.
                
                % Decode DATA field
                PSDU = this.PLCP_DATA_Decode(PPDU_cfo_cor, h_est, start_of_payload_samples, rxvector, noise_variance, SNR_est_avg);
            end
            

        end 
        
        
        
        
    end
    
    methods (Access = protected)
        
       
       function preamble = PLCP_Preamble_Generate(this)
            
           % Call superclass method first
            preamble = PLCP_Preamble_Generate@PLCP_OFDM(this); 
            
            %Remember: PLCP_Preamble_Generate() in PLCP_OFDM class repeats
            % N_Tx times LTS. So for MIMO transmissions (where N_Tx=2),
            % there is 2 LTS copy concatenated in the preamble.
            % Delete the corresponding LTS copy from the preamble to give the
            % necessary orthogonality to the MIMO transmission of LTS.
            % Which LTS copy to delete depends on stream tx_id.
            preamble(160*(this.sts_repeat_factor+mod(this.tx_id,this.SysVar.N_Tx)) : 160*(this.sts_repeat_factor+1+mod(this.tx_id,this.SysVar.N_Tx))-1) = zeros(1,160);
            
        end
        
       function [decision, start_of_payload] = PLCP_LTS_Signal_Detection(this, preamble)
            %Frame/signal detection using LTS. In MIMO transmissions, there is
            %more than one LTS field, depending on the number of
            %transmitters. This function considers this extension to 802.11a preamble 
            %for the final decision of start of payload, and detection.
            
            end_of_lts = zeros(1,this.SysVar.N_Tx);
            
            %For each LTS field, cancel the most powerful symbol and
            %re-correlate to find more peaks.
            for i = 1:this.SysVar.N_Tx
            
                %Complex cross correlation of Rx waveform with time-domain LTS
                lts_corr = abs(conv(conj(fliplr(this.lts_t)), sign(preamble)));
                
                %There will be ideally 2*N_Tx peaks after cross correlation with the
                %original LTS symbol.
                
                
                %Find all correlation peaks
                lts_peaks = find(lts_corr > this.SysVar.LTSCorr_THRESH*max(lts_corr));
                
                %Select best candidate correlation peak as LTS-payload boundary
                [LTS1, LTS2] = meshgrid(lts_peaks,lts_peaks);
                [lts_second_peak_indices,y] = find(LTS2-LTS1 == length(this.lts_t));
                
                %No valid correlation peak found
                if(isempty(lts_second_peak_indices) || lts_peaks(max(lts_second_peak_indices)) < 161)
                    break;
                else
                    
                    %Set the sample indices of the payload symbols and preamble
                    end_of_lts(i)  = lts_peaks(max(lts_second_peak_indices));
                
                    preamble(end_of_lts(i)-160:end_of_lts(i)-1)=zeros(1,160);
                end
                

            end
            
            
            if(~isempty(end_of_lts))
                
                if(this.SysVar.N_Tx == 2)
                    % In case of N_Tx==2, we check if we found the peaks from
                    % the second L-LTF field or not?
                    % If it is from second, it should be somewhere after the
                    % STSs + 1st LTS samples (ideally). Plus we leave 80 smpls of fudge factor.
                    %
                    % If not from the second L-LTF symbol, then add 160 samples
                    % more to find the start of payload.
                    
                    if(max(end_of_lts) > this.sts_repeat_factor+1*160+80)
                        start_of_payload = max(end_of_lts);
                    else
                        start_of_payload = max(end_of_lts)+160;
                         warning('Could not find the second the LTS field rom cross correlation. \n');
                    end
                    
                    decision=true;
                else
                    error('PLCP_OFDM_MIMO->PLCP_LTS_Signal_Detection: N_Tx != 2 is not supported yet!');
                end
            else
                decision = false;
                start_of_payload = -1;
                this.ErrStats.add_DetectionErr();
                fprintf('LTS Detection Failed!\n');

            end
       end
        
       function out = PLCP_Coarse_CFO (this, raw_rx, start_of_payload)
           % Coarse (Time domain) CFO estimation and correction.
   
           if(this.NLTS_Samples ~= this.SysVar.N_Tx*160)
               error('PLCP_OFDM_MIMO->PLCP_Coarse_CFO:There is smth wrong about the size of LTS samples generated!');
           end
           
           for i=1:this.SysVar.N_Tx
               
               lts_ind = start_of_payload - i*160; % There is 160 samples of L-LTF for each transmitter.
               
               %Extract LTS (not yet CFO corrected)
               rx_lts  = raw_rx(lts_ind : lts_ind+159);
               rx_lts1 = rx_lts(-64+-this.SysVar.FFT_OFFSET + [97:160]);
               rx_lts2 = rx_lts(-this.SysVar.FFT_OFFSET + [97:160]);
               
               %Estimate coarse CFO
               rx_cfo_est_lts = median(unwrap(angle(rx_lts1 .* conj(rx_lts2))));
               rx_cfo_est_lts = rx_cfo_est_lts/(2*pi*64);
               cfo_est_lts(mod(i,2)+1) =  rx_cfo_est_lts;
               %fprintf('[MIMO] LTS - (RX #%d) CFO Est: %3.2f kHz\n',mod(i,2)+1, rx_cfo_est_lts*20e6*1e-3);
           end
           
           %Correction 
           % - average of CFOs from all LTSs for each TX.
           rx_cfo_corr_t = exp(1i*2*pi*mean(cfo_est_lts)*[0:length(raw_rx)-1]');
           out = raw_rx .* rx_cfo_corr_t;
        
           %Test for STS time domain estimation for checking.
           %rx_sts = raw_rx(1:160); rx_sts1 = rx_sts(1:16); rx_sts2 = rx_sts(17:32); rx_sts3 = rx_sts(33:48); rx_sts4 = rx_sts(49:64); cfo_sts = mean(unwrap(angle(rx_sts3.*conj(rx_sts4))))/(2*pi*16);cfo_sts*20e6*1e-3
           
       end
       
       function data_cfo_corr = PLCP_FREQ_CFO_MIMO_Estimate (this, data_f, h_est, NOFDM_sym)
           % Frequency domain CFO estimation for MIMO signals (multiple
           % output-multiple input). This assumes distributed pilots among
           % tx streams in frequency). 
           
           %Initialize local variables
           eq_p_ant= zeros(length(this.index_subc_pilots),NOFDM_sym);
           data_symbols = data_f(this.index_subc_data,:,:);
            
             %Equalize pilots per each antenna (Remember pilots are shared
             %among TX streams. We also take it into account). 
             for rx=1:this.Nr
                 
                 for tx=1:this.SysVar.N_Tx
                     eq_p_ant(tx:this.SysVar.N_Tx:end,:) = data_f(this.index_subc_pilots(tx:this.SysVar.N_Tx:end),:,rx)./ repmat(h_est(this.index_subc_pilots(tx:this.SysVar.N_Tx:end),tx,rx),1,NOFDM_sym);
                 end
                 
                 %Find the average phase error per antenna from pilots
                 pilots_f_mat = repmat(this.pilots_f_all.', 1, NOFDM_sym);
                 pilot_phase_err_ant= angle(mean(eq_p_ant.*pilots_f_mat));
                 
                 pilot_cfo_estimate_ant = unwrap(pilot_phase_err_ant);
                 
                 pilot_phase_corr_ant = repmat(exp(-1i*pilot_cfo_estimate_ant), this.NOFDM_data_subc, 1);
                 
                 %Apply the pilot phase correction per symbol
                 data_cfo_corr(:,:,rx) = data_symbols(:,:,rx) .* pilot_phase_corr_ant;
             end
                 
           
       end
       
       function PSDU = PLCP_DATA_Decode(this, rx_samples_t_cor, h_est, start_of_payload, rxvector, noise_var, snr)
             %Received DATA field samples processing
    
             %Extract necessary parameters from RXVECTOR
            n_dbps = this.Ndbps(1,this.RateID(rxvector.DATARATE)); % Number of data bits per symbol
            n_cbps = this.Ncbps(1,this.RateID(rxvector.DATARATE)); % Number of coded bits per symbol
            n_bpsc = this.Nbpsc(1,this.RateID(rxvector.DATARATE)); % Number of coded bits per subcarrier
            CR = this.CodingRate(1,this.RateID(rxvector.DATARATE)) ; % Code Rate
            M=this.ConstellationSize(1,this.RateID(rxvector.DATARATE)); % Constellation Size
            hdemod_soft = modem.qamdemod('M',M,'SymbolOrder','Gray','OutputType','Bit','DecisionType','approximate llr'); % Create Modulator object
            hdemod = modem.qamdemod('M',M,'SymbolOrder','Gray','OutputType','Bit','DecisionType','hard decision'); % Create Modulator object
            const_scale  = modnorm(hdemod.Constellation,'avpow',1); %  Scale the constellation for average power equal to 1.
          

            
            %Calculate the total number of OFDM symbols & padding bits
            NOFDM_sym = ceil((16 + 8*rxvector.LENGTH + 6)/n_dbps);
            Npad = NOFDM_sym * n_dbps - (16 + 8*rxvector.LENGTH + 6);
            
           
             %Extract the DATA samples (# OFDM symbol+GI)
             start_of_DATA = start_of_payload(:) + (this.NOFDM_FFT_size+this.CPLength); %One OFDM symbol skip due to (SIGNAL)
             end_of_DATA = start_of_DATA(:)+NOFDM_sym*(this.NOFDM_FFT_size+this.CPLength)-1;
             
             if(any(end_of_DATA>size(rx_samples_t_cor,1)))
                 fprintf('PLCP_OFDM_MIMO->PLCP_DATA_Decode: Not enough samples for DATA field. \n Probably wrong SIGNAL-LENGTH field or DATARATE estimation!\n');
                 this.ErrStats.add_SIGNALNVErr();
                 PSDU=-1;
                 return;
             end
             
             if(this.SysVar.short_GI==0)
                 for i=1:this.Nr
                    data_samples(:,i) = rx_samples_t_cor(start_of_DATA(i) : end_of_DATA(i),i);
                 end
                 
                 data_mat = reshape(data_samples, (this.NOFDM_FFT_size+this.CPLength), NOFDM_sym,[]);
                 
                 %Remove the cyclic prefix, keeping FFT_OFFSET samples of CP (on average)
                 data_noCP = data_mat(this.CPLength-this.SysVar.FFT_OFFSET+[1:this.NOFDM_FFT_size], :,:);
             else
                 error('PLCP_OFDM_MIMO->PLCP_DATA_Decode: Short_GI is not supported!');
             end
             
             %Take the FFT
             data_f = fft(data_noCP, this.NOFDM_FFT_size);

             h_est_datasc = h_est(this.index_subc_data,:,:); %Only data subcarriers.
                          
    
             %Frequency domain CFO estimation (for distributed MIMO tx)           
             data_symbols_cfo_cor = this.PLCP_FREQ_CFO_MIMO_Estimate(data_f, h_est, NOFDM_sym);
             
             % MIMO Decoding
             
             %Check if Maximum Likelihood (ML) Decoder is chosen for
             %decoding. Can be slow for 64-QAM modulations.
             if(this.SysVar.ML_decoder)
                 data_bits_scrambled = this.PLCP_MIMO_ML_Decoder(data_symbols_cfo_cor, h_est_datasc, hdemod_soft, NOFDM_sym,n_cbps,n_bpsc, CR, noise_var, snr);
             else
                 data_bits_scrambled = this.PLCP_MIMO_Decoder(data_symbols_cfo_cor, h_est_datasc, hdemod_soft, NOFDM_sym,n_cbps,n_bpsc, CR, noise_var, snr);
             end
             
             for tx_stream=1:this.SysVar.N_Tx
                     
                 try
                     %De-scramble
                     data = step(this.Descrambler,data_bits_scrambled(:,tx_stream));
                     %Construct SERVICE field and PSDU
                     % Optional: Check if SERVICE is all zero?if so, set RXVECTOR.
                     %            Remove padding and tail bits
                     
                     service = data(1:16);
                     PSDU(:,tx_stream) = data(17:end-6-Npad);
                     
                 catch ME
                     fprintf('PLCP_OFDM_MIMO->PLCP_DATA_Decode: Size of scrambled input : %d\n', size(data_bits_scrambled(:,tx_stream),1));
                     this.ErrStats.add_SIGNALNVErr();
                     PSDU=-1;
                     return;
                 end

             end
   
          
       end
        
       function data_bits_scrambled = PLCP_MIMO_ML_Decoder(this, data_symbols_cfo_cor, h_est_datasc, hdemod_soft, NOFDM_sym,n_cbps,n_bpsc, CR, noise_var, snr)
            %ML(Maximum Likelihood) MIMO decoder
            
            hmod =  modem.qammod(hdemod_soft); %Modulator    
            modOrd = log2(hmod.M); %Modulation order
            const_scale  = modnorm(hdemod_soft.Constellation,'avpow',1); %  Scale the constellation for average power equal to 1.
            
            allBits  = de2bi(0:2^(modOrd*this.Nt)-1, 'left-msb')';
            allTxSig = const_scale*reshape(modulate(hmod, allBits(:)), this.Nt, 2^(modOrd*this.Nt));
            
            %Initialize variables
            llr_all = zeros(this.NOFDM_data_subc,modOrd*NOFDM_sym,this.Nt);
            
            for subc = 1:1:this.NOFDM_data_subc
                H = squeeze(h_est_datasc(subc,:,:)).';
                llr_per_subc = [];
                for ofdm = 1:NOFDM_sym
                    r = squeeze(data_symbols_cfo_cor(subc,ofdm,:));
                    [lambdaML, k] = min(sum(abs(repmat(r,[1,2^(modOrd*this.Nt)]) - H*allTxSig).^2));
                    symML = allTxSig(:,k); %ML estimated symbol combination
                    bitsML = allBits(:,k);
                    
                    %Compute the log-likelihood ratios (LLR) for each bit
                    for jb_integer=1:length(bitsML)
                        complement_symbols = allTxSig(:,find(allBits(jb_integer,:)==~bitsML(jb_integer,1))); 
                        lambdaML_complement = min(sum(abs(repmat(r,[1,size(complement_symbols,2)]) - H*complement_symbols).^2));
                        
                        if(bitsML(jb_integer)==0)
                            llr(jb_integer,1) = lambdaML-lambdaML_complement;
                        else
                            llr(jb_integer,1) = lambdaML_complement - lambdaML;
                        end
                    end
                    
                    llr_per_subc = vertcat(llr_per_subc,reshape(llr,[],this.Nt));
                    
                end
                
                llr_all(subc,:,:) = llr_per_subc;
            end
            
            
            
            for next_stream=1:this.Nt
                data_bits_soft=[];
                data_bits_datasc_soft=llr_all(:,:,next_stream);
            % Reshape llr_all so that output bits are in
            % order as: same symbol bits->subcarrier->ofdm symbol
                for i=1:NOFDM_sym
                    f = data_bits_datasc_soft(:,(i-1)*log2(hmod.M)+1:i*log2(hmod.M));
                    v= reshape(f.',[],1);
                    data_bits_soft=vertcat(data_bits_soft,v);
                end
                
            %Strangely MATLAB computes LLR as probability of bit being 0 versus bit being 1
            % Demodulator soft output and decoder soft-input assumes the
            % same
            data_bits_soft = -1*data_bits_soft; 
            
            % De-interleave the soft bits (LLRs)
            % - Nbpsc=log2(M),Ncbps=this.NOFDM_data_subc
            data_deintlv = deinterleaver(data_bits_soft,16,n_cbps,n_bpsc);
            
            % Decode the bits (using soft decision viterbi decoding)
            % Traceback table length = 3; Not specified in standards.
            
            %For different code rates: 
            try
                if CR == 1/2
                    % No puncturing
                    data_bits_scrambled(:,next_stream) = vitdec(data_deintlv,this.Trellis,3, 'trunc','unquant');
                    
                elseif CR == 2/3
                    % Puncturing pattern [1 1 1 0] for 2/3 rate
                    data_bits_scrambled(:,next_stream) = vitdec(data_deintlv,this.Trellis,3, 'trunc','unquant',[1 1 1 0]);
                    
                elseif CR == 3/4
                    % Puncturing pattern [1 1 1 0 0 1] for 3/4 rate
                    data_bits_scrambled(:,next_stream) = vitdec(data_deintlv,this.Trellis,3, 'trunc','unquant',[1 1 1 0 0 1]);
                    
                else
                    fprintf('PLCP_OFDM_MIMO->PLCP_MIMO_ML_Decoder: Code Rate(CR) argument is not valid.\n');
                    this.ErrStats.add_SIGNALNVErr();
                    PSDU=-1;
                    return;
                end
            catch ME
                data_deintlv
                error(ME.message);
            end
                
                
            end
       end
       
       function data_bits_scrambled = PLCP_MIMO_Decoder(this, data_symbols_cfo_cor, h_est_datasc, hdemod_soft, NOFDM_sym,n_cbps,n_bpsc, CR, noise_var, snr)
            %MIMO decoder with/without successive interference cancellation.
            % MMSE-SIC, ZF-SIC
            
            %Initialize local variables
            est_stream_sym = zeros(this.NOFDM_data_subc,NOFDM_sym,this.Nr);
            total_streams_decoded = 0; %Total number of spatial streams processed.
            
            const_scale  = modnorm(hdemod_soft.Constellation,'avpow',1); %  Scale the constellation for average power equal to 1.
            M = hdemod_soft.M; % Modulation size
            hmod =  modem.qammod(hdemod_soft); %Modulator
            
            % Optimal Stream Ordering based on SNR
            [~,order] = sort(mean(snr),'descend'); %Avg SNR for each TX stream
            
            % Start decoding streams using the order specified.
            for next_stream = order
                
                data_bits_soft=[];%Empty vector to collect ordered soft data bits.
                total_streams_decoded = total_streams_decoded + 1;
                
                if(this.SysVar.SIC_decoder)
                    %If SIC decoder is used, remove the previous stream
                    %from rest of the signal first.
                    data_symbols_cfo_cor = data_symbols_cfo_cor - est_stream_sym;
                    
                    %If it is the last stream to be decoded, we can use MRC
                    if(total_streams_decoded==length(order))
                       
                        %Initialize numerator and denominator structures for MRC
                        %combining.
                        Num = zeros(this.NOFDM_data_subc,NOFDM_sym);
                        Denom = zeros(this.NOFDM_data_subc,NOFDM_sym);
                        
                        for i = 1:this.Nr
                            % Maximum Ratio combining of received signal
                            %
                            % x_hat = (conj(h)*y / (h^H*h)) = (conj(h)*y / sum(abs(h)^2))
                            %
                            Num = Num + conj(repmat(squeeze(h_est_datasc(:,next_stream,i)),1,NOFDM_sym)).* data_symbols_cfo_cor(:,:,i);
                            Denom = Denom + abs(repmat(squeeze(h_est_datasc(:,next_stream,i)),1,NOFDM_sym)).^2;
                            
                        end
                        
                        x_decorr =  Num./Denom;
                    end
                end
                
                for subc = 1:this.NOFDM_data_subc
                    
                    if((this.SysVar.SIC_decoder && total_streams_decoded < length(order)) || (~this.SysVar.SIC_decoder))
                        %Decorrelator
                        y = squeeze(data_symbols_cfo_cor(subc,:,:)).';
                        w = pinv(squeeze(h_est_datasc(subc,:,:)).');
                        for ofdm = 1:NOFDM_sym
                            x_decorr(subc,ofdm) = w(next_stream,:)*y(:,ofdm);
                        end
                    end
                    
                    %Demodulate symbols (Output soft outputs (LLR))
                    data_bits_datasc_soft(subc,:) = reshape(demodulate(hdemod_soft,(1/const_scale)*reshape(x_decorr(subc,:),[],1)),1,[]);
                    
                end
                
                
            % Reshape data_bits_datasc_soft so that output bits are in
            % order as: same symbol bits->subcarrier->ofdm symbol
            for i=1:NOFDM_sym
                f = data_bits_datasc_soft(:,(i-1)*log2(M)+1:i*log2(M));
                v= reshape(f.',[],1);
                data_bits_soft=vertcat(data_bits_soft,v);
            end
            
            % De-interleave
            % - Nbpsc=log2(M),Ncbps=this.NOFDM_data_subc
            data_deintlv = deinterleaver(data_bits_soft,16,n_cbps,n_bpsc);
            
            % DECODER (Soft decision decoding)
            % Viterbi Decoder
            % Traceback table length = 3; Not specified in standards.
            try
                if CR == 1/2
                    % No puncturing
                    data_bits_scrambled(:,next_stream) = vitdec(data_deintlv,this.Trellis,3, 'trunc','unquant');
                    
                elseif CR == 2/3
                    % Puncturing pattern [1 1 1 0] for 2/3 rate
                    data_bits_scrambled(:,next_stream) = vitdec(data_deintlv,this.Trellis,3, 'trunc','unquant',[1 1 1 0]);
                    
                elseif CR == 3/4
                    % Puncturing pattern [1 1 1 0 0 1] for 3/4 rate
                    data_bits_scrambled(:,next_stream) = vitdec(data_deintlv,this.Trellis,3, 'trunc','unquant',[1 1 1 0 0 1]);
                    
                else
                    fprintf('PLCP_OFDM_MIMO->PLCP_MIMO_Decoder: Code Rate(CR) argument is not valid.\n');
                    this.ErrStats.add_SIGNALNVErr();
                    PSDU=-1;
                    return;
                end
            catch ME
                data_deintlv
                error(ME.message);
            end
            
            % Successive Interference Cancellation (SIC)
            % ------------------------------------------------------
            % (If SIC_decoder==1, we try cancel to cancel decoded streams
            % from the rest of the signal.)
            % NOTE: SIC can propogate errors which significantly effect the
            % decoder performance. If decoding time is not an issue, we
            % suggest to use ML Decoder implementation in this class. 
            
            
            if(this.SysVar.SIC_decoder && (total_streams_decoded <length(order)))
                       
                % Encode the scrambled data (Section 18.3.5.6 - Page 1597-1598)
                % Encoding the DATA field for re-creating the signal to be
                % cancelled
                data_coded=convenc(data_bits_scrambled(:,next_stream),this.Trellis); 
                % Puncturing (--Need to puncture if the rate is not 1/2)
                switch CR
                    case 2/3
                        data_coded(4:4:end)=[]; % Remove every 4th bit
                    case 3/4
                        data_coded(4:6:end)=[];
                        data_coded(4:5:end)=[]; % Remove every 4th+5th bit
                end
                
                % Data interleaving (Section 18.3.5.7 -- Page 1598)
                %   (Ncbps=1(BPSK),Nbpsc=this.NOFDM_data_subc)
                data_intlv = interleaver(data_coded,16,n_cbps,n_bpsc);
                
                % Modulate 
                tx_symbols = modulate(hmod,data_intlv);
                tx_symbols = const_scale*tx_symbols; %All in column vector.
                
                % The stream of complex numbers is divided into groups
                % of NOFDM_data_subc = 48 complex numbers.
                tx_symbols_matrix(1:this.NOFDM_data_subc,1:NOFDM_sym) = reshape(tx_symbols,this.NOFDM_data_subc,NOFDM_sym);
                
                
                %Add the effects of the channel in frequency domain 
                % This is necessary for the proper cancellation of the
                % stream from the rest of the signal.
                for rx_ant = 1:this.Nr
                    est_stream_sym(:,:,rx_ant) = repmat(h_est_datasc(:,next_stream,rx_ant),1,NOFDM_sym).*tx_symbols_matrix;
                end
            end
            
            
            end
            
        end
    end
    
end

