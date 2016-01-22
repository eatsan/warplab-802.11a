%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% © All rights reserved. ECOLE POLYTECHNIQUE FEDERALE DE LAUSANNE, Switzerland, 
% Laboratory Algorithmic Research in Network Information (ARNI), 2015
% AUTHORS: Emre ATSAN
%
% See the LICENSE.TXT file for more details.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


classdef PLCP_OFDM < PLCP
    %PLCP_OFDM OFDM PHY PLCP implementation class
    %
    %   Inherits PLCP abstract class.
    %
    % INFO: This class implements the SISO PLCP for OFDM PHY in IEEE
    % 802.11a OFDM PHY standard (see IEEE 802.11-2012 standard document). 
    
    
    properties
        SysVar                  % System variables struct
        RateID                  % Transmission rates mapping (RateValue -> ID)
        RateTable               % Transmission rates table (Rate Value)
        Rate4Bit                % Rate values to TXVECTOR bit representation
        CodingRate              % Available code rates vector
        ConstellationSize       % Available const. sizes vector (BPSK, 4-QAM, 16-QAM, 64-QAM)
        NOFDM_FFT_size=64;      % Number of OFDM subcarriers/FFT size
        NOFDM_data_subc=48;     % Number of OFDM data subcarriers
        NOFDM_pilot_subc=4;     % Number of OFDM pilot subcarriers
        Nbpsc                   % SET / Number of coded bits per subcarrier
        Ncbps                   % SET / Number of coded bits per OFDM symbol
        Ndbps                   % SET / Number of data bits per OFDM symbol
        index_subc_data         % Indices of OFDM data subcarriers
        index_subc_pilots       % Indices of OFDM pilot subcarriers
        index_subc_null         % Indices of OFDM null subcarriers
        index_subc_used         % Indices of OFDM subcarriers which are used
        CPLength                % Cyclic Prefix Length
        Trellis                 % Encoder/decoder trellis
        ModV_BPSK               % Fast Modulation Vector for BPSK (only in SIGNAL block)
        ModF_BPSK               % Fast Modulation Function for BPSK
        DeModF_BPSK             % Fast Demodulation Function for BPSK
        Scrambler               % DATA scrambler
        Descrambler             % DATA descrambler
        NPreamble_Samples       % Number of Preamble Samples 
        NLTS_Samples            % Number of LTS field Samples 
        lts_f                   % Long training sequnce (in frequency domain).
        lts_t                   % Long training symbols (in time domain)
        sts_t                   % Short training symbols (in time domain)
        lts_t_ext               % Extended LTS (in time domain)
        sts_t_ext               % Extended STS (in time domain) (to fit 8musec frame)
        sts_repeat_factor       % Number of sts_t_ext OFDM symbol repetition for better AGC adjustment
        pilots_f                % Pilot frequency domain sequence for this node. 
                                % -- (Used to support distributed pilot transmissions in Distributed-MIMO implementation)
        pilots_f_all            % Complete Pilot frequency domain sequence.
        ErrStats                % PLCP layer Error statistics object
        Nt                      % Number of TX chains available
        Nr                      % Number of RX chains available
        
        %For debugging
        CoarseCFOEst
        SNR         % SNR estimation for each stream per antenna
        start_of_payload % Where payload starts after synchronization.
    end
    
    % Public methods (accessible outside the class)
    methods (Access = public)
        
        % Constructor
        function p=PLCP_OFDM(sysvar)
            
            p.SysVar = sysvar;
            p.Nt = p.SysVar.N_Tx;
            p.Nr = p.SysVar.N_Rx;
            
            % Set class parameters according to channel bandwidth used
            if(sysvar.Ch_BW==20)
                % Table 18-6, Page 1595
                p.RateID = containers.Map([6, 9, 12, 18, 24, 36, 48, 54],1:8);
                p.RateTable = [6, 9, 12, 18, 24, 36, 48, 54];
                p.Rate4Bit = [[1,1,0,1];[1,1,1,1];[0,1,0,1];[0,1,1,1];[1,0,0,1];[1,0,1,1];[0,0,0,1];[0,0,1,1]];
                
            elseif (sysvar.Ch_BW==10) 
                % TODO: Support for 10 MHz channels are not implemented. 
                % Table 18-6, Page 1595
                p.RateID = containers.Map([3, 4.5, 6, 9, 12, 18, 24, 27],1:8);
                p.RateTable = [6, 9, 12, 18, 24, 36, 48, 54];
                p.Rate4Bit = [[1,1,0,1];[1,1,1,1];[0,1,0,1];[0,1,1,1];[1,0,0,1];[1,0,1,1];[0,0,0,1];[0,0,1,1]];
            else
                error('PLCP_OFDM: Unsupported Channel Bandwidth for the constructor.')
            end
            
            % Set class parameters according to short/long guard interval
            % (GI) options. 
            if(sysvar.short_GI==0)
                p.CPLength = p.NOFDM_FFT_size/4; % 0.8 msec
            elseif (sysvar.short_GI==1)
                p.CPLength = p.NOFDM_FFT_size/8; % 0.4 msec
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Modulation-dependent parameters
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Table 18-4, Page 1590
            p.CodingRate = [1/2,3/4,1/2,3/4,1/2,3/4,2/3,3/4];
            p.ConstellationSize = [2,2,4,4,16,16,64,64]; % BPSK=2, QPSK=4, 16-QAM=16, 64-QAM=64
            p.Nbpsc = [1,1,2,2,4,4,6,6]; %Number of coded bits per subcarrier
            p.Ncbps = p.NOFDM_data_subc.*p.Nbpsc; % Number of coded bits per OFDM symbol
            p.Ndbps =  p.CodingRate.*p.Ncbps; %Number of data bits per OFDM symbol
            
            % Convolutional Code (Constraint length 7, generator poly [133, 171]
            % Page 1597 - Section 18.3.5.6, IEEE 802.11-2012.
            p.Trellis = poly2trellis(7,[133 171]); % Encoder/decoder trellis
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Define OFDM subcarrier mask for identification of data, pilots, and null
            % subcarriers
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % This mask is definied following the 802.11-2012 standard pgs
            % 1603-1604.
            % I use 0 to label subcarriers that are set tozero, 1 to label suncarriers
            % that carry payload, and 2 to label subcarriers tha carry data
            % The ordering of the vector entries is [subcarrier -32, subcarrier -31, ...., subcarrier -1, subcarrier 0, subcarrier 1, ...., subcarrier 30, subcarrier 31];
            % ORIGINAL AS PER 802.11:
            subc_ID_mask =  [0 0 0 0 0 0 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2 0 2 2 2 2 2 2 1 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 0 0 0 0 0];
            
            % Apply fft shift
            subc_ID_mask = fftshift(subc_ID_mask);
            
            % Get indices of data and pilot subcarriers
            p.index_subc_data = find(subc_ID_mask == 2);
            p.index_subc_pilots = find(subc_ID_mask == 1);
            p.index_subc_null = find(subc_ID_mask == 0);
            p.index_subc_used = union(p.index_subc_data, p.index_subc_pilots);
            
            
            %Fast Modulation Vector & Function for BPSK
            p.ModV_BPSK  =  [-1 1];
            p.ModF_BPSK = @(x) complex(p.ModV_BPSK(1+x),0);
            p.DeModF_BPSK = @(x) double(real(x)>0);
            
            %Create the scrambler & descrambler to be used for DATA block
            % -- based on the 802.11-2012 standard.
            p.Scrambler = comm.Scrambler(2,[0 -4 -7],[1 0 1 1 1 0 1]);
            p.Descrambler = comm.Descrambler(2,[0 -4 -7],[1 0 1 1 1 0 1]);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Preamble training samples
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Short training symbols (STS)
            p.sts_repeat_factor = 2;
            sts_f = zeros(1,64);
            sts_f(1:27) = [0 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0];
            sts_f(39:64) = [0 0 1+1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0];
            p.sts_t = ifft(sqrt(13/6).*sts_f, 64); %Time-domain STS.
            p.sts_t_ext = repmat([p.sts_t p.sts_t p.sts_t(1:32)],1,p.sts_repeat_factor); % Extension of STS to fit 8musec OFDM frame.
            
            
            %Long training symbols (LTS)
            % -- for CFO and channel estimation
            p.lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
            p.lts_t = ifft(p.lts_f, 64); %Time domain LTS signal;
            p.lts_t_ext = [p.lts_t(33:64) p.lts_t p.lts_t]; % Extension of the LTS symbols
            
            % Number of LTS samples 
            %  Note that if there is MIMO transmission (e.g.N_Tx>1), 
            %  then there is an orthogonal LTS training per tx antenna.
            p.NLTS_Samples = sysvar.N_Tx*length(p.lts_t_ext);
            
            %Number of preamble samples: LTS + STS samples
            p.NPreamble_Samples = length(p.sts_t_ext)+p.NLTS_Samples ;
            
            % Pilot subcarrier values (TODO: implement pps)
            p.pilots_f_all = [1,1,-1,1];
            p.pilots_f = [1,1,-1,1]; %Single transmitter case, all pilots assigned to node.
            
            %Create an PLCP_ErrorStats object to be used in this class.
            p.ErrStats = PLCP_ErrorStats();
            
        end
        
        
        % PLCP RX procedure
        function [rxvector,PSDU] = PLCP_RX(this, PPDU)
            %Inputs:
            % PPDU: PLCP Protocol Data Unit. Input from PMD. 
            %
            %Output:
            % PSDU: PLCP Service Data Unit (MAC payload - MPDU).
            % rxvector: RXVECTOR parameters for MAC layer operations.
            
            % Check if the number of received antennas matches with the
            % received raw signal sets from the lower layers (PMD).
            if(this.Nr ~= size(PPDU,2))
                warning('PLCP_OFDM->PLCP_RX: Not enough sample dimensions received from the PMD layer!');
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
            
            
            if(all(decision)==false)
                warning('PLCP_OFDM_SIMO->PLCP_RX: Frame Detection failed on all antennas!');
                rxvector=-1; PSDU=-1;
                return;
            end
            
            
            rx_stream = 1;
            
            for rx_antenna=find(decision)
                
                %This variable was used to check if it is better to decide
                %using the stronger antenna value for all antennas. However
                %it did not work out, so just ignore it for now. 
                master_antenna =rx_antenna;
                
                %Sanity Check for synchronization
                if ((start_of_payload_samples(master_antenna)-this.NLTS_Samples) < -5) %More than 5 samples off is not a good estimate.
                    warning(['PLCP_OFDM_SIMO->PLCP_RX: Frame Detection on receive antenna (',rx_antenna,') failed! Wrong start of payload estimation.']);
                    decision(rx_antenna) = false;
                    continue;
                elseif ((start_of_payload_samples(master_antenna)-this.NLTS_Samples) < 1) %Less than 5 samples off can be corrected.
                    start_of_payload_samples(master_antenna) = start_of_payload_samples(master_antenna) + 1 - (start_of_payload_samples(master_antenna)-this.NLTS_Samples);
                end
                
                % Coarse CFO estimation & correction
                PPDU_cfo_cor(:,rx_stream) = this.PLCP_Coarse_CFO(squeeze(PPDU(:,rx_antenna)),start_of_payload_samples(master_antenna), rx_antenna); %time domain
                
                %Channel Estimation (using LTS training)
                preamble_lts = PPDU_cfo_cor((start_of_payload_samples(master_antenna)-this.NLTS_Samples):(start_of_payload_samples(master_antenna)-1),rx_antenna);
                h_est(:,:,rx_stream) = this.PLCP_Channel_Estimate(preamble_lts); % Frequency domain channel estimates
                
                
                % SNR estimate (using LTS training)
                [SNR_est(:,:,rx_stream),SNR_est_avg(rx_stream,:), noise_variance(:,:,rx_stream) ]= this.PLCP_SNR_est_LTS(preamble_lts);
                
                %Set the start of payload samples only among the detected
                %antennas
                start_p_samples_detected(rx_stream) = start_of_payload_samples(rx_antenna);
                
                rx_stream = rx_stream + 1; %Move to the next stream to process.
                
            end
            
            
% DEBUGGING / Channel Estimates            
%             if(this.Nr==2)
%                 %Plot channel magnitudes (with fftshift)
%                 h_shifted1 = fftshift(h_est(:,1,1));
%                 h_shifted2 = fftshift(h_est(:,1,2));
%                 figure; plot(0.3125:0.3125:20,abs(h_shifted1(:)),'-*'); hold on; plot(0.3125:0.3125:20,abs(h_shifted2(:)),'-*r');grid on;
%                 title('Channel magnitude per subcarrier (with FFTshift');
%                 xlabel('Baseband Bandwidth (MHz)');
%                 legend('RXA','RXB');
%                 drawnow;
%             end
            
             this.SNR = SNR_est_avg; 
             this.start_of_payload = start_of_payload_samples;
             
            % Extract SIGNAL field & Create RXVECTOR
            [length, d_rate] = this.PLCP_SIGNAL_Decode(PPDU_cfo_cor, h_est, start_p_samples_detected);
            
            if(length==-1 && d_rate==-1)
                warning('PLCP_OFDM->PLCP_RX: SIGNAL not decoded!');
                rxvector=-1; PSDU=-1;
                return;
            else
                rxvector = RXVECTOR_OFDMPHY(length, d_rate);
                % TODO: Ideally this needs to create an event to signal MAC Layer.
                %  -- PHY-RXSTART.indication(RXVECTOR) primitive.
                
                % Decode the DATA field from received samples
                PSDU = this.PLCP_DATA_Decode(PPDU_cfo_cor, h_est, start_p_samples_detected, rxvector, noise_variance);
            end
            
            
        end
        
        % PLCP TX procedure
        function PPDU = PLCP_TX(this, txvector, PSDU)
            %Inputs:
            % txvector: TXVECTOR parameters for PHY operations.
            % PSDU: PLCP Service Data Unit (MAC payload).
            %
            %Output:
            % PPDU: PLCP Protocol Data Unit. Input to PMD.
            
            %Generate preamble
            preamble = this.PLCP_Preamble_Generate();
            
            %Generate PLCP SIGNAL field
            signal = this.PLCP_SIGNAL_Generate(txvector.LENGTH, txvector.DATARATE);
            
            %Generate PLCP DATA field
            data = this.PLCP_DATA_Generate(txvector, PSDU);
            
            %Create the PPDU 
            PPDU = [preamble, signal, data];
        end
        
        
    end
    
    methods (Access = protected)
        
        function preamble = PLCP_Preamble_Generate(this)
            
            
            % Use 10 copies of the STS + 2.5 copies of LTS
            % (see Page 1593, Std. 802.11-2012)
            if(this.SysVar.Ch_BW==20) % Only implemented for 20Mhz channels with standard GI2.
                if(this.SysVar.short_GI==0)
                    %For MIMO case, use N_Tx LTS transmission 
                    %(for orthogonal channel estimation)
                    preamble = [this.sts_t_ext  repmat(this.lts_t_ext, 1, this.SysVar.N_Tx)];
                else
                    error('PLCP_OFDM->PLCP_Preamble_Generate: Short_GI is not supported!');
                end
            else
                error('PLCP_OFDM->PLCP_Preamble_Generate: Only 20Mhz channel bandwidth supported.');
            end
            
        end
        
        function samples_t = PLCP_SIGNAL_Generate(this, length, rate)
            %SIGNAL OFDM symbol will be encoded with a Conv. Encoder of
            %rate 1/2 and bits are modulated with BPSK.
            
            signal_f=zeros(1,24); %24 bit for SIGNAL field (Page 1595)
            
            
            if(this.SysVar.Ch_BW==20)
                signal_f(1:4) = this.Rate4Bit(this.RateID(rate),:); %rate field
                signal_f(5) = 0; %Reserved bit for future use.
                signal_f(6:17) =  fliplr(sscanf(dec2bin(length/8,12),'%1d').'); %length field (should be in bytes).
                signal_f(18) = mod(sum(signal_f(1:17)),2); % even parity bit
                % signal_f(19:24) = zeros(1,6); SIGNAL TRAIL is already set to zero.
                
                %Encode the signal field with conv encoder
                signal_coded=convenc(signal_f,this.Trellis); % Encoding the SIGNAL field
                
                %Interleave the bits (Ncbps=1(BPSK),Nbpsc=this.NOFDM_data_subc)
                signal_intlv = interleaver(signal_coded,16,this.NOFDM_data_subc,1);
                
                %Modulate the coded signal field bits with BPSK
                symbols_f(this.index_subc_data) = arrayfun(this.ModF_BPSK, signal_intlv);
                symbols_f(this.index_subc_pilots) = this.pilots_f; %TODO: implement polarity sequence
                
                %OFDM symbol generation
                samples_t = ifft(symbols_f, this.NOFDM_FFT_size);
                
                % Adding Guard interval
                if(this.SysVar.short_GI==0)
                    samples_t = appendGI(samples_t,this.CPLength); %appendGI returns back row vector.
                    samples_t = reshape(samples_t,1,[]); % Reshape back to column vector.
                else
                    error('PLCP_OFDM->PLCP_SIGNAL_Generate: Short_GI is not supported!');
                end
            else
                error('PLCP_OFDM->PLCP_SIGNAL_Generate: Only 20Mhz channel options are supported.');
            end
        end
        
        function samples_t = PLCP_DATA_Generate(this, txvector, PSDU)
            %DATA field is generated according to specs (Section 18.3.5, Page 1596)
            
            %Sanity Check
            if(txvector.LENGTH ~= length(PSDU))
                error('Verify your MAC layer PSDU size and TXVECTOR.LENGTH field matches!')
            else
                psdu_length = txvector.LENGTH; % NOTE: in bits.
            end
            
            % Extract necessary parameters from TXVECTOR
            n_dbps = this.Ndbps(1,this.RateID(txvector.DATARATE)); % Number of data bits per symbol
            n_cbps = this.Ncbps(1,this.RateID(txvector.DATARATE)); % Number of coded bits per symbol
            n_bpsc = this.Nbpsc(1,this.RateID(txvector.DATARATE)); % Number of coded bits per subcarrier
            CR = this.CodingRate(1,this.RateID(txvector.DATARATE)) ; % Code Rate
            M=this.ConstellationSize(1,this.RateID(txvector.DATARATE)); % Constellation Size
            hmod = modem.qammod('M',M,'SymbolOrder','Gray','InputType','Bit'); % Create Modulator object
            const_scale  = modnorm(hmod.Constellation,'avpow',1); %  Scale the constellation for average power equal to 1.
            
            
            % Extend the resulting bit string with zero bits (at least 6 bits)
            % so that the resulting length is a multiple of NDBPS. The
            % resulting bit string constitutes the DATA part of the packet.
            % Refer to 18.3.5.4 for details.
            data_len=ceil((psdu_length+16+6)/n_dbps)*n_dbps;
            
            % 16bits from SERVICE, 6 bits from PPDU TAIL, rest is padding.
            data_f = zeros(data_len,1); %Initialize DATA field.
            
            %Number of OFDM symbols 
            NOFDM_sym = ceil((16 + psdu_length + 6)/n_dbps);
            %Number of padding bits
            Npad = NOFDM_sym * n_dbps - (16 + psdu_length + 6);
            
            %fprintf('--  DATA: # of Pad bits = %4.0f --- # of OFDM sym: %4.0f\n',Npad, NOFDM_sym);
            
            % Append the PSDU to the SERVICE field of the TXVECTOR.
            data_f(1:16,1) = txvector.SERVICE;
            data_f(17:16+psdu_length,1) = PSDU;
            
            % Scramble the data (Section 18.3.5.5 - Page 1597)
            % (Note: Pointless if you have a random data bit generation for test).
            data_scrambled = step(this.Scrambler,data_f);
            
            % Encode the scrambled data (Section 18.3.5.6 - Page 1597-1598)
            data_coded=convenc(data_scrambled,this.Trellis); % Encoding the DATA field
            
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
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Consider possible PMD Layer relocation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Modulate the DATA field bits with correspoing BPSK or QAM
            % constellation.
            tx_symbols = modulate(hmod,data_intlv);
            tx_symbols = const_scale*tx_symbols; %All in column vector.
            
            % The stream of complex numbers is divided into groups
            % of NOFDM_data_subc = 48 complex numbers.
            tx_symbols_matrix(1:this.NOFDM_data_subc,1:NOFDM_sym) = reshape(tx_symbols,this.NOFDM_data_subc,NOFDM_sym);
            
            symbols_f_matrix = zeros(this.NOFDM_FFT_size,NOFDM_sym); %Initialize the matrix input to OFDM modulator.
            
            symbols_f_matrix(this.index_subc_data,1:NOFDM_sym) = tx_symbols_matrix; % Assign to data subcarriers
            symbols_f_matrix(this.index_subc_pilots,1:NOFDM_sym) = repmat(this.pilots_f.',1,NOFDM_sym); %TODO: implement polarity sequence
            
            %OFDM symbol generation
            symbols_t_matrix = ifft(symbols_f_matrix, this.NOFDM_FFT_size);
            samples_t_matrix = zeros(this.NOFDM_FFT_size+this.CPLength, NOFDM_sym);
            
            % Adding Guard interval
            if(this.SysVar.short_GI==0)
                for sym=1:NOFDM_sym %For each OFDM symbol
                    samples_t_matrix(:,sym) = appendGI(symbols_t_matrix(:,sym),this.CPLength);
                end
                samples_t = reshape(samples_t_matrix,1,[]);
            else
                error('PLCP_OFDM->PLCP_DATA_Generate: Short_GI is not supported!');
            end
        end
        
        function out = PLCP_Coarse_CFO (this, raw_rx, start_of_payload, rx_antenna)
            % Coarse (Time domain) CFO estimation and correction.
            
            if(this.NLTS_Samples ~= 160)
                error('There is smth wrong about the size of LTS samples generated!');
            end
            
            lts_ind = start_of_payload - 160; % There is 160 samples of LTS for SISO.

            % Extract LTS (Coarse CFO corrected)
            rx_lts  = raw_rx(lts_ind : lts_ind+159);
            %rx_lts  = raw_rx_coarse_corr(lts_ind : lts_ind+159);
            rx_lts1 = rx_lts(-64+-this.SysVar.FFT_OFFSET + [97:160]);
            rx_lts2 = rx_lts(-this.SysVar.FFT_OFFSET + [97:160]);
            
            % Coarse CFO (with mean estimate)
            rx_fine_cfo_est_lts = median(unwrap(angle(rx_lts1 .* conj(rx_lts2))));
            rx_fine_cfo_est_lts = rx_fine_cfo_est_lts/(2*pi*64);
            cfo_lts_value = -rx_fine_cfo_est_lts*20e6*1e-3;
            
            %fprintf('[RX %d] LTS CFO Est : %3.4f kHz\n', rx_antenna, cfo_lts_value);
            
            %Correction
            rx_cfo_corr_t = exp(1i*2*pi*rx_fine_cfo_est_lts*[0:length(raw_rx)-1]');
            out = raw_rx .* rx_cfo_corr_t;
       


%  For DEBUGGING    
%             
%             if(rx_antenna==1)
%                 this.CoarseCFOEst(1)=cfo_lts_value;
%                 figure;
%                 plot((-1*20e6*1e-3/(2*pi*64))*unwrap(angle(rx_lts1 .* conj(rx_lts2))),'-*b');
%             else
%                 this.CoarseCFOEst(2)=cfo_lts_value;
%                 figure;
%                 plot((-1*20e6*1e-3/(2*pi*64))*unwrap(angle(rx_lts1 .* conj(rx_lts2))),'-*r');
%             end
%             grid on;
%             xlabel('Sample No');
%             ylabel('CFO estimate (kHz)');
%             drawnow

%            rx_lts_corr = out(lts_ind : lts_ind+159);
%            rx_lts1_corr = rx_lts_corr(-64+-this.SysVar.FFT_OFFSET + [97:160]);
%            rx_lts2_corr = rx_lts_corr(-this.SysVar.FFT_OFFSET + [97:160]);
%            rx_cfo_est_lts_corr = mean(unwrap(angle(rx_lts1_corr .* conj(rx_lts2_corr))));
%            rx_cfo_est_lts_corr = rx_cfo_est_lts_corr/(2*pi*64);
%            fprintf('LTS CFO Est (after correction): %3.2f kHz\n', -rx_cfo_est_lts_corr*20e6*1e-3);
        end
        
        function [decision,start_of_data] = PLCP_Signal_Detection(this, preamble)
            % Signal Detection from the STS/LTS (training seq).
            %
            % The first six symbols of the short preamble are for:
            % signal detection, automatic gain control(AGC), and
            % diversity selection (for antenna).
            % The last four short symbols are intended
            % for coarse frequency offset estimation (Time CFO) and timing
            % synchronization.
            %
            %             %Complex cross correlation of Rx waveform with time-domain STS
            %             % - We use first 16 bit of the STS OFDM symbol.
            %             sts_corr = abs(conv(conj(fliplr(this.sts_t(1:16))), sign(preamble)));
            %
            %             %Skip early and late samples
            %             %sts_corr = sts_corr(15:end-160);
            %
            %             % Find all correlation peaks
            %             sts_peaks = find(sts_corr > this.SysVar.LTSCorr_THRESH*max(sts_corr));
            %
            %
            %             if(length(sts_peaks)==this.sts_repeat_factor*10)
            %                 decision=true;
            %                 start_of_data = sts_peaks(this.sts_repeat_factor*10)+this.SysVar.N_Tx*length(this.lts_t_ext)+1;
            %             else
            %                 decision=false;
            %                 start_of_data = -1;
            %                 fprintf('STS Detection Failed!\n');
            %                 this.ErrStats.add_STSDetectionErr();
            %             end
            %
            
            % IMPORTANT NOTE: I have disabled synchronization from STS, 
            % since AGC takes some time in WARP radios. 
            % Also since we have all the frame samples
            % already received, STS detection is not speeding anything
            % anyways. If you want to try STS detection, uncomment the code
            % above and delete the next line.
            decision=false;
            
            %If the STS signal detection failed, try LTS detection:
            if(decision==false)
                % cross-correlation based frame
                % detection using LTS signals. This part is optional,
                % increases load at the receiver, however gives better
                % detection performance.
                [decision, start_of_data]  = this.PLCP_LTS_Signal_Detection(preamble);
                
            end
            
        end
        
        function [decision, start_of_payload] = PLCP_LTS_Signal_Detection(this, preamble)
            % Frame detection using the LTS field of the signal. This
            % technique is not a standard way of signal detection, however
            % we try this if Frame detection failed using STS only.
            %
            % This function only works for *single* L-LTF field transmission
            % in the preamble. If there is more than one L-LTF concatenated
            % in the preamble, use PLCP_OFDM_MIMO subclass corresponding
            % function for detection.
            
            %Complex cross correlation of Rx waveform with time-domain LTS
            lts_corr = abs(conv(conj(fliplr(this.lts_t)), sign(preamble)));
            
            %Skip early and late samples
            lts_corr = lts_corr(32:end-32);
            
            %Find all correlation peaks
            lts_peaks = find(lts_corr > this.SysVar.LTSCorr_THRESH*max(lts_corr));
            
            %Select best candidate correlation peak as LTS-payload boundary
            [LTS1, LTS2] = meshgrid(lts_peaks,lts_peaks);
            [lts_second_peak_index,y] = find(LTS2-LTS1 == length(this.lts_t));
            
            %No valid correlation peak found
            if(isempty(lts_second_peak_index))
                fprintf('LTS Detection Failed!\n');
                
                decision=false;
                start_of_payload = -1;
                this.ErrStats.add_DetectionErr();
                return;
            end
            
            %Set the sample indices of the payload symbols and preamble
            start_of_payload  = lts_peaks(max(lts_second_peak_index))+32;
            decision=true;
            
        end
        
        function h_est = PLCP_Channel_Estimate(this, preamble_lts)
            
            %Initialize the channel matrix estimate 
            h_est = zeros(this.NOFDM_FFT_size,this.SysVar.N_Tx);
            
            for i=1:this.SysVar.N_Tx
                
                % Start of LTS preamble index for TX antenna mod(i,this.SysVar.N_Tx)+1.
                %  --Observe that we first estimate the second LTS field
                
                % In case of MIMO tx, there is 160 samples of L-LTF 
                % for each transmitter.
                lts_ind = length(preamble_lts)+1 - i*160; 
                
                rx_lts = preamble_lts(lts_ind : lts_ind+159);
                %Re-extract LTS for channel estimate
                rx_lts1 = rx_lts(-64+-this.SysVar.FFT_OFFSET + [97:160]);
                rx_lts2 = rx_lts(-this.SysVar.FFT_OFFSET + [97:160]);
                
                rx_lts1_f = fft(rx_lts1); %LTS-1 in freq domain
                rx_lts2_f = fft(rx_lts2); %LTS-2 in freq domain
                
                %Calculate channel estimate (averaging over time for
                %degrading noise effects)
                h_est(:,mod(i,this.SysVar.N_Tx)+1) = this.lts_f.' .* ((rx_lts1_f + rx_lts2_f)/2);
                
            end
            
        end
        
        function [ SNR_hat_2_dB, SNR_hat_2_avg_dB, W_hat_per_subc ] = PLCP_SNR_est_LTS( this, preamble_lts )
            %PLCP_SNR_est_LTS Estimate Estimate channel SNR from the LTS
            
            %Set to 1 to get SNR information printed to command window for
            % every transmission.
            DEBUG =0 ;
           
            
            % Generate error if it is not a 160 sample per Tx node preamble lts.
            assert(length(preamble_lts)==this.SysVar.N_Tx*160);
            
            for i=1:this.SysVar.N_Tx
                
                % Start of LTS preamble index for TX antenna mod(i,this.SysVar.N_Tx)+1.
                %  --Observe that we first estimate the second LTS field
                
                % In case of MIMO tx, there is 160 samples of L-LTF 
                % for each transmitter.
                lts_ind = length(preamble_lts)+1 - i*160; 
                
                y = preamble_lts(lts_ind : lts_ind+159);
                
                % Extract time domain lts symbols (2 64 samples after GI removed).
                % lts1=y(33:32+64); lts2 = y(32+64+1:32+64+64);
                lts1 = y(-64+-this.SysVar.FFT_OFFSET + [97:160]);
                lts2 = y(-this.SysVar.FFT_OFFSET + [97:160]);
                
                
                % FFT
                lts1_f_all = fft(lts1);lts2_f_all = fft(lts2);
                
                % SNR estimation method is based on the following paper:
                % (SNR Estimation Algorithm Based on the preamble 
                %   for OFDM Systems in Frequency Selective Channels, G. Ren et al.)
                y_1_all = lts1_f_all .* conj(this.lts_f.');
                y_2_all = lts2_f_all .* conj(this.lts_f.');
                
                
                % Only choose used/data subcarriers (Ideally all the carriers, but somehow in
                % sim, it always see 0 for unused subcarriers).
                y_1 = y_1_all(this.index_subc_data); y_2 = y_2_all(this.index_subc_data);
                % Noise variance estimate per subcarrier
                W_hat_per_subc(:,mod(i,this.SysVar.N_Tx)+1) = 1/2 * abs(y_1 - y_2).^2;
                
                % Noise variance estimate average over data subcarriers
                W_hat(mod(i,this.SysVar.N_Tx)+1) = mean(W_hat_per_subc(:,mod(i,this.SysVar.N_Tx)+1));
                
                %Per subcarrier SNR
                SNR_hat_2 = (((abs(y_1).^2 + abs(y_2).^2)/ 2 ) / W_hat(mod(i,this.SysVar.N_Tx)+1) ) - 1 ; 
                %Per subcarrier SNR in dB
                SNR_hat_2_dB(:,mod(i,this.SysVar.N_Tx)+1) = 10*log10(SNR_hat_2); 
                % Average SNR over used subcarriers
                SNR_hat_2_avg = mean(SNR_hat_2); 
                % Average SNR over used subcarriers in dB
                SNR_hat_2_avg_dB(mod(i,this.SysVar.N_Tx)+1) = 10*log10(SNR_hat_2_avg);
                
                if(DEBUG)
                    fprintf('[Spatial # %d] Avg. SNR estimate from LTS: %3f dB.\n',mod(i,this.SysVar.N_Tx)+1, SNR_hat_2_avg_dB(mod(i,this.SysVar.N_Tx)+1));
                end
                
            end

        end
        
        function [data_length, data_rate] = PLCP_SIGNAL_Decode(this, rx_samples_t_cor, h_est, start_of_payload)
            % This function decodes the SIGNAL field in a PLCP frame
            
            Num = zeros(this.NOFDM_FFT_size,1);
            Denom = zeros(this.NOFDM_FFT_size,1);
            active_rx_chains = size(h_est,3);
            
            if(active_rx_chains~=size(rx_samples_t_cor,2))
                error('PLCP_OFDM->PLCP_SIGNAL_Decode: Channel Estimate dimensions and received vector dimensions do not match!')
            end
            
            for rx_antenna=1:active_rx_chains
                
                %Extract the SIGNAL samples (One OFDM symbol+GI)
                if(this.SysVar.short_GI==0)
                    signal_samples = rx_samples_t_cor(start_of_payload(rx_antenna) : start_of_payload(rx_antenna)+(this.NOFDM_FFT_size+this.CPLength)-1,rx_antenna);
                    signal_mat = reshape(signal_samples, (this.NOFDM_FFT_size+this.CPLength), 1);
                    
                    %Remove the cyclic prefix, keeping FFT_OFFSET samples of CP (on average)
                    signal_noCP = signal_mat(this.CPLength-this.SysVar.FFT_OFFSET+[1:this.NOFDM_FFT_size], :);
                else
                    error('PLCP_OFDM->PLCP_SIGNAL_Decode: Short_GI is not supported!');
                end
                %Take the FFT
                signal_f = fft(signal_noCP, this.NOFDM_FFT_size);
                
                %During multiple TX antenna (in single TX it is h_est)
                %
                %Sum up the total channel effects since both transmitters
                %send the same SIGNAL field. (TODO: Investigate the validity
                %of this)
                h_est_combined = sum(h_est(:,:,rx_antenna),2);
                
                % SIGNAL PILOTS in dist. mimo transmission are orthogonal,
                % so better to estimate the channel w/o summing up the 
                % independent channel estimates for only shared pilots.
                for tx_ant=1:this.Nt
                    h_est_combined(this.index_subc_pilots(tx_ant:2:end)) = h_est(this.index_subc_pilots(tx_ant:2:end),tx_ant,rx_antenna);
                end
                
                
                %Frequency domain CFO estimation for the signal field
                signal_f_corr = this.PLCP_FREQ_CFO_Estimate(signal_f(:,:),h_est_combined, 1);
                
                
                % Will be used for maximum ratio combining
                Num = Num + conj(h_est_combined).* signal_f_corr(1:1:this.NOFDM_FFT_size,1);
                Denom = Denom + abs(h_est_combined).^2;
                
            end
            
            
            %Channel Equalization (maximum ratio combining - MRC)
            signal_eq_f = Num(1:1:this.NOFDM_FFT_size,1)./Denom(1:1:this.NOFDM_FFT_size,1);
            
            
            signal_symbols_datasc = signal_eq_f(this.index_subc_data,:);
            
            % DEBUG signal symbols
            %figure;
            %scatterplot(signal_symbols_datasc); title('SIGNAL field symbols');grid on;
            
            %Demodulate symbols (SIGNAL field is always BPSK modulated)
            signal_bits_datasc = arrayfun(this.DeModF_BPSK, signal_symbols_datasc);
            
            %De-interleave the bits
            % - Ncbps=1(BPSK),Nbpsc=this.NOFDM_data_subc
            signal_deintlv = deinterleaver(signal_bits_datasc,16,this.NOFDM_data_subc,1);
            
            %Decode (Hard decision decoding)
            % TODO: Soft decision decoding
            % Traceback table length = 3; Not specified in standards.
            signal_bits = vitdec(signal_deintlv,this.Trellis,3, 'trunc','hard');
            
            
            %Check parity for integrity
            % If parity passed, extract necessary fields for RXVECTOR.
            if(signal_bits(18) ~= mod(sum(signal_bits(1:17)),2)) % even parity bit
                fprintf('SIGNAL field parity bit failed!\n');
                this.ErrStats.add_SIGNALNVErr();
                data_length=-1;
                data_rate=-1;
                return;
            end
            
            %Extract the data length field
            data_length = bin2dec(num2str(fliplr(signal_bits(6:17).')));
            
            %Assertion
            if(data_length < 1 && data_length > this.SysVar.maxLength)
                fprint('SIGNAL Length field decoding failed!\n');
                this.ErrStats.add_SIGNALNVErr();
                data_length=-1;
                data_rate=-1;
                return;
            end
            
            %Extract the data rate bits
            rate_bits = signal_bits(1:4).';
            %Find the corresponding rate id from data rate bits.
            rate_id = find(ismember(this.Rate4Bit,rate_bits,'rows'),1);
           
            if(isempty(rate_id))
                fprintf('SIGNAL Datarate bits are not meaningful!\n');
                this.ErrStats.add_SIGNALNVErr();
                data_length=-1;
                data_rate=-1;
                return;
            else
                %Set the tx rate for the DATA field 
                data_rate = this.RateTable(rate_id);
            end
            
            tail = signal_bits(19:24);%should be six zeros
            if (sum(tail)~=0)%should be zero
                fprintf('SIGNAL Tail bits are non-zero!\n');
                this.ErrStats.add_SIGNALNVErr();
                data_length=-1;
                data_rate=-1;
                return;
            end
            
            
            
        end
        
        function PSDU = PLCP_DATA_Decode(this, rx_samples_t_cor, h_est, start_of_payload, rxvector, noise_var)
            %Received DATA field samples processing
            
            %Extract necessary parameters from RXVECTOR
            n_dbps = this.Ndbps(1,this.RateID(rxvector.DATARATE)); % Number of data bits per symbol
            n_cbps = this.Ncbps(1,this.RateID(rxvector.DATARATE)); % Number of coded bits per symbol
            n_bpsc = this.Nbpsc(1,this.RateID(rxvector.DATARATE)); % Number of coded bits per subcarrier
            CR = this.CodingRate(1,this.RateID(rxvector.DATARATE)) ; % Code Rate
            M=this.ConstellationSize(1,this.RateID(rxvector.DATARATE)); % Constellation Size
            
            % Create Modulator object
            hdemod = modem.qamdemod('M',M,'SymbolOrder','Gray','OutputType','Bit','DecisionType','approximate llr');
            const_scale  = modnorm(hdemod.Constellation,'avpow',1); %  Scale the constellation for average power equal to 1.
            
            
            %Calculate the total number of OFDM symbols & padding bits
            NOFDM_sym = ceil((16 + 8*rxvector.LENGTH + 6)/n_dbps);
            Npad = NOFDM_sym * n_dbps - (16 + 8*rxvector.LENGTH + 6);
            
            %Extract the DATA samples (# OFDM symbol+GI)
            start_of_DATA = start_of_payload + (this.NOFDM_FFT_size+this.CPLength); %One OFDM symbol skip due to (SIGNAL)
            end_of_DATA = start_of_DATA+NOFDM_sym*(this.NOFDM_FFT_size+this.CPLength)-1;
            
            %Find how many RX chains is active from received structure
            active_rx_chains = size(h_est,3);
            
            %Initialize numerator and denominator structures for MRC
            %combining.
            Num = zeros(this.NOFDM_FFT_size,NOFDM_sym);
            Denom = zeros(this.NOFDM_FFT_size,NOFDM_sym);
            
            %Initialize noise variance matrix for MRC transformed signal
            noise_var_per_subc= zeros(this.NOFDM_data_subc,1);
            
            %Empty vector to collect ordered soft data bits.
            data_bits_soft=[];
            
            % Sanity check for signal dimension
            if(active_rx_chains~=size(rx_samples_t_cor,2))
                error('PLCP_OFDM->PLCP_DATA_Decode: Channel Estimate dimensions and received vector dimensions do not match!')
            end
            
            % Sanity check for received signal size and END_OF_DATA
            % estimation.
            if(all(end_of_DATA>length(rx_samples_t_cor)))
                fprintf('PLCP_OFDM->PLCP_DATA_Decode: Not enough samples for DATA field. \n Probably wrong SIGNAL-LENGTH field or DATARATE estimation!\n');
                this.ErrStats.add_SIGNALNVErr();
                PSDU=-1;
                return;
            end
            
            
            if(this.SysVar.short_GI==0)
                
                %Extract DATA field samples
                for i=1:active_rx_chains
                    data_samples(:,i) = rx_samples_t_cor(start_of_DATA(i) : end_of_DATA(i),i);
                end
                
                data_mat = reshape(data_samples, (this.NOFDM_FFT_size+this.CPLength), NOFDM_sym,[]);
                
                %Remove the cyclic prefix, keeping FFT_OFFSET samples of CP (on average)
                data_noCP = data_mat(this.CPLength-this.SysVar.FFT_OFFSET+[1:this.NOFDM_FFT_size], :,:);
            else
                error('PLCP_OFDM->PLCP_DATA_Decode: Short_GI is not supported!');
            end
            
            % Take the FFT
            data_f = fft(data_noCP, this.NOFDM_FFT_size);
            
            
            % MRC (Maximum Ratio Combining) in case of multiple receive
            % chains
            for i=1:active_rx_chains
                
                %Frequency domain CFO estimation
                % -> Observe that we do CFO estimation before channel equalization!
                %  --> In this case, we need to do ch. equalization for
                %  pilots inside CFO_Estimate function.
                data_f_corr = this.PLCP_FREQ_CFO_Estimate(data_f(:,:,i),h_est(:,1,i),NOFDM_sym);
                
                % Maximum Ratio combining of received signal
                %
                % x_hat = (conj(h)*y / (h^H*h)) = (conj(h)*y / sum(abs(h)^2))  
                %
                Num = Num + conj(repmat(squeeze(h_est(:,:,i)),1,NOFDM_sym)).* data_f_corr(:,:);
                Denom = Denom + abs(repmat(squeeze(h_est(:,:,i)),1,NOFDM_sym)).^2;
                
                % Noise variance transformation to MRC noise.
                %
                % Since we are doing MRC for the received signals from
                % different antennas, the final combined signal has a
                % different noise characteristic than the originals. 
                
                % MRC basically combines the received signals with a
                % weighted sum, where the weights are the channel
                % (h_est(i)) powers for each antenna i. Similarly, we can
                % take the weighted sum of the noise variances per antenna
                % for the MRC combined noise_variance that will be used
                % with the soft demodulator.
                
                % Exact transformation is as follows: (W_i is the noise
                % variance computed from RX chain i).
                % W_mrc = W_1/abs(h1)^2 + W_2/abs(h2)^2
                
                noise_var_per_subc = noise_var_per_subc + (noise_var(:,:,i))./ abs(squeeze(h_est(this.index_subc_data,:,i))).^2;
            end
            
            
            %MRC Equalization
            data_eq_f = Num(1:1:this.NOFDM_FFT_size,:)./Denom(1:1:this.NOFDM_FFT_size,:);
            data_symbols_datasc = data_eq_f(this.index_subc_data,:);
            
            
            % DEBUG: For plotting the received data symbols
            if(this.SysVar.plot_data_symbols)
             symbols_plot_update(reshape(data_symbols_datasc,[],1), hdemod);
            end
            
            % Per subcarrier based demodulation.
            % The idea is that some subcarriers have a higher noise
            % variance (especially due to sampling clock harmonics at ch. 6 & 14. 
            % Using average noise variance instead of per subcarrier can 
            % affect LLR computation and eventually the soft decoder performance.
            % NOTE: This is very very experimental. 
            for subc=1:this.NOFDM_data_subc
                hdemod.NoiseVariance = (1/const_scale)*noise_var_per_subc(subc,1);
                %Demodulate symbols (Output soft outputs (LLR))
                data_bits_datasc_soft(subc,:) = reshape(demodulate(hdemod,(1/const_scale)*reshape(data_symbols_datasc(subc,:),[],1)),1,[]);
            end
            
            % Reshape data_bits_datasc_soft so that output bits are in
            % order: bits from same symbol -> subcarrier -> ofdm symbol
            for i=1:NOFDM_sym
                f = data_bits_datasc_soft(:,(i-1)*log2(M)+1:i*log2(M));
                v= reshape(f.',[],1);
                data_bits_soft=vertcat(data_bits_soft,v);
            end
            
            % De-interleaver (as in IEEE 802.11-2012 standard.)
            % - Nbpsc=log2(M),Ncbps=this.NOFDM_data_subc
            data_deintlv = deinterleaver(data_bits_soft,16,n_cbps,n_bpsc);
            
            % DECODER (Soft decision decoding)
            % Viterbi Decoder
            % Traceback table length = 3; Not specified in standards.
            try
                if CR == 1/2
                    % No puncturing
                    data_bits = vitdec(data_deintlv,this.Trellis,3, 'trunc','unquant');
                    
                elseif CR == 2/3
                    % Puncturing pattern [1 1 1 0] for 2/3 rate
                    data_bits = vitdec(data_deintlv,this.Trellis,3, 'trunc','unquant',[1 1 1 0]);
                    
                elseif CR == 3/4
                    % Puncturing pattern [1 1 1 0 0 1] for 3/4 rate
                    data_bits = vitdec(data_deintlv,this.Trellis,3, 'trunc','unquant',[1 1 1 0 0 1]);
                    
                else
                    fprintf('PLCP_OFDM->PLCP_DATA_Decode: Code Rate in RXVECTOR is not valid.\n');
                    this.ErrStats.add_SIGNALNVErr();
                    PSDU=-1;
                    return;
                end
            catch ME
                data_deintlv
                error(ME.message);
            end
            
            % De-scrambler
            data_f = step(this.Descrambler,data_bits);
            
            %Construct SERVICE field and PSDU
            % Optional: Check if SERVICE is all zero? if so, set RXVECTOR.
            %            Remove padding and tail bits
            
            service = data_f(1:16);
            PSDU = data_f(17:end-6-Npad); % TAIL & Padding removed
            
            
        end
        
        function signal_data_cfo_cor = PLCP_FREQ_CFO_Estimate(this, signal_f, h_est, nOFDM_sym)
            %Frequency domain CFO estimation from pilot subcarriers.
            %
            % Inputs: This function can handle several OFDM symbols if the fft
            % output is reshaped into [# of subcarriers x # of OFDM sym ]
            %
            % Returns: Freq. domain CFO corrected data  (matrix)
            
            %Repeat the pilots across all OFDM symbols
            pilots_f_mat = repmat(this.pilots_f.', 1, nOFDM_sym);
            
            % Extract the pilots and calculate per-symbol phase error
            % Since this function should be run before channel
            % equalization, we do a simple ZF equalization for the pilots.
            pilots_rx_f = signal_f(this.index_subc_pilots, :)./repmat(h_est(this.index_subc_pilots),1,nOFDM_sym);
            pilot_phase_err_initial = angle(pilots_rx_f.*pilots_f_mat);
            
            numeratorweight = abs(h_est(this.index_subc_pilots));
            numeratorweight = reshape(numeratorweight,[],1);
            weight_matrix = zeros(length(this.index_subc_pilots),nOFDM_sym);
            weight_matrix(:,:) = repmat(numeratorweight,1,nOFDM_sym);
            
            
            % Average the per pilot phase difference using the weighted average
            freqEst = sum(weight_matrix .* pilot_phase_err_initial) ./ sum(weight_matrix);
            
            pilot_phase_corr = repmat(exp(-1i*freqEst), this.NOFDM_FFT_size, 1);
            
            %Apply the pilot phase correction per symbol
            signal_data_cfo_cor = signal_f .* pilot_phase_corr;
            
            
            
        end
        
    end
    
end