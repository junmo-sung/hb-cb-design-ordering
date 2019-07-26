function [ SP ] = main_function( SP )
SP.Mt = SP.Nt/SP.Lt * SP.sym;
SP.Mr = SP.Nr/SP.Lr;
SP.M  = SP.Mt*SP.Mr;
% SP.Gr = SP.G; SP.Gt = SP.G;
% SP.NtBeam           = SP.Nt;
% SP.NrBeam           = SP.Nr;
SP.NtBlck           = SP.Nt / SP.Lt;
SP.NrBlck           = SP.Nr / SP.Lr;
SP.SNR_lin_array    = 10.^(.1*SP.SNR_db_array);
% Channel-related parameters
SP.virtualChannel   = 'spatial';
SP.integerTapDelay  = 0; % Is each tap delay an integer?
SP.BeamAlign        = 0; % AoA and AoD aligned with dictionary grid
SP.Nc               = 1; % # of channel taps
SP.beta             = 0.8; % a roll-off factor of the shaping filter

% NMSE result matrices initialization
SP.NMSE_LS_sung     = zeros(SP.iter, length(SP.SNR_db_array));
SP.NMSE_LS_lee      = zeros(SP.iter, length(SP.SNR_db_array));
SP.NMSE_LS_md       = zeros(SP.iter, length(SP.SNR_db_array));
SP.NMSE_OMP_sung    = zeros(SP.iter, length(SP.SNR_db_array));
SP.NMSE_OMP_lee     = zeros(SP.iter, length(SP.SNR_db_array));
SP.NMSE_OMP_md      = zeros(SP.iter, length(SP.SNR_db_array));
SP.NMSE_BPDN_sung   = zeros(SP.iter, length(SP.SNR_db_array));
SP.NMSE_BPDN_lee    = zeros(SP.iter, length(SP.SNR_db_array));
SP.NMSE_BPDN_md     = zeros(SP.iter, length(SP.SNR_db_array));
SP.NMSE_IHT_sung    = zeros(SP.iter, length(SP.SNR_db_array));
SP.NMSE_BPDN_lee    = zeros(SP.iter, length(SP.SNR_db_array));
SP.NMSE_IHT_md      = zeros(SP.iter, length(SP.SNR_db_array));
SP.NMSE_OR          = zeros(SP.iter, length(SP.SNR_db_array));
SP.SE_opt           = zeros(SP.iter, length(SP.SNR_db_array));
SP.SE_oracle        = zeros(SP.iter, length(SP.SNR_db_array));
SP.SE_LS_sung       = zeros(SP.iter, length(SP.SNR_db_array));
SP.SE_LS_lee        = zeros(SP.iter, length(SP.SNR_db_array));
SP.SE_LS_md         = zeros(SP.iter, length(SP.SNR_db_array));
SP.SE_OMP_sung      = zeros(SP.iter, length(SP.SNR_db_array));
SP.SE_OMP_lee       = zeros(SP.iter, length(SP.SNR_db_array));
SP.SE_OMP_md        = zeros(SP.iter, length(SP.SNR_db_array));
SP.SE_BPDN_sung     = zeros(SP.iter, length(SP.SNR_db_array));
SP.SE_BPDN_lee      = zeros(SP.iter, length(SP.SNR_db_array));
SP.SE_BPDN_md       = zeros(SP.iter, length(SP.SNR_db_array));
SP.SE_IHT_sung      = zeros(SP.iter, length(SP.SNR_db_array));
SP.SE_IHT_lee       = zeros(SP.iter, length(SP.SNR_db_array));
SP.SE_IHT_md        = zeros(SP.iter, length(SP.SNR_db_array));

% Beamformer generation

[BF.sung.Phibar,    BF.sung.PhibarNorm, BF.sung.Sigma,  BF.sung.Phi,    BF.sung.AtBar,  BF.sung.ArBar,  BF.sung.WBarDiag,  BF.sung.gridBeamPairs] = sung(SP);

% SNR iteration
for snrInd = 1:length(SP.SNR_db_array)
    SP.rho = SP.SNR_lin_array(snrInd)*SP.sigma^2;
    % Simulation iteration
    for simInd = 1:SP.iter
        fprintf('Simulation INFO:  SNR=%d (%d/%d), Iteration=%d/%d \n', SP.SNR_db_array(snrInd), snrInd, length(SP.SNR_db_array), simInd, SP.iter)
        % Beamformer generation
        [BF.md.Phibar,      BF.md.PhibarNorm,   BF.md.Sigma,    BF.md.Phi,      BF.md.AtBar,    BF.md.ArBar,    BF.md.Q_diag,      BF.md.gridBeamPairs  ] = random_cb(SP);
        %         if strcmp(SP.rfArchitecture, 'PS')
        [BF.lee.Phibar,     BF.lee.PhibarNorm,  BF.lee.Sigma,   BF.lee.Phi,     BF.lee.AtBar,   BF.lee.ArBar,   BF.lee.WBarDiag,   BF.lee.gridBeamPairs] = mtc_cb(SP);
        %         end
        
        % Channel instance
        [SP.H, SP.A_t, SP.A_r, SP.beamPairs] = channelGen(SP);
        
        % Noise instance
        SP.noiseVector = SP.sigma*sqrt(1/2) * (randn(SP.M*SP.Nr, 1) + 1j*randn(SP.M*SP.Nr, 1));
        
        % Call algorithms
        [Hest_OMP_sung,  Hest_BPDN_sung,  Hest_IHT_sung,  Hest_LS_sung,   Hest_oracle,    SP.NMSE_LS_sung(simInd, snrInd),   SP.NMSE_OMP_sung(simInd, snrInd), SP.NMSE_BPDN_sung(simInd, snrInd),   SP.NMSE_IHT_sung(simInd, snrInd),    SP.NMSE_OR(simInd, snrInd)] = main_sung(SP, BF.sung);
        [Hest_OMP_lee,   Hest_BPDN_lee,   Hest_IHT_lee,   Hest_LS_lee,    ~,              SP.NMSE_LS_lee(simInd, snrInd),    SP.NMSE_OMP_lee(simInd, snrInd),  SP.NMSE_BPDN_lee(simInd, snrInd),    SP.NMSE_IHT_lee(simInd, snrInd),     SP.NMSE_OR(simInd, snrInd)] = main_mtc(SP, BF.lee);
        [Hest_OMP_md,    Hest_BPDN_md,    Hest_IHT_md,    Hest_LS_md,                     SP.NMSE_LS_md(simInd, snrInd),     SP.NMSE_OMP_md(simInd, snrInd),   SP.NMSE_BPDN_md(simInd, snrInd),     SP.NMSE_IHT_md(simInd, snrInd),                                ] = main_random(SP, BF.md);
        
        % Spectral Efficiency
        SP.SE_opt           (simInd, snrInd) = SE( SP.H, SP.H,              SP, SP.SNR_lin_array(snrInd) );
        SP.SE_oracle        (simInd, snrInd) = SE( SP.H, Hest_oracle,       SP, SP.SNR_lin_array(snrInd) );
        SP.SE_LS_sung       (simInd, snrInd) = SE( SP.H, Hest_LS_sung,      SP, SP.SNR_lin_array(snrInd) );
        SP.SE_LS_md         (simInd, snrInd) = SE( SP.H, Hest_LS_md,        SP, SP.SNR_lin_array(snrInd) );
        SP.SE_OMP_sung      (simInd, snrInd) = SE( SP.H, Hest_OMP_sung,     SP, SP.SNR_lin_array(snrInd) );
        SP.SE_OMP_md        (simInd, snrInd) = SE( SP.H, Hest_OMP_md,       SP, SP.SNR_lin_array(snrInd) );
        SP.SE_BPDN_sung     (simInd, snrInd) = SE( SP.H, Hest_BPDN_sung,    SP, SP.SNR_lin_array(snrInd) );
        SP.SE_BPDN_md       (simInd, snrInd) = SE( SP.H, Hest_BPDN_md,      SP, SP.SNR_lin_array(snrInd) );
        SP.SE_IHT_sung      (simInd, snrInd) = SE( SP.H, Hest_IHT_sung,     SP, SP.SNR_lin_array(snrInd) );
        SP.SE_IHT_md        (simInd, snrInd) = SE( SP.H, Hest_IHT_md,       SP, SP.SNR_lin_array(snrInd) );
        SP.SE_LS_lee        (simInd, snrInd) = SE( SP.H, Hest_LS_lee,       SP, SP.SNR_lin_array(snrInd) );
        SP.SE_OMP_lee       (simInd, snrInd) = SE( SP.H, Hest_OMP_lee,      SP, SP.SNR_lin_array(snrInd) );
        SP.SE_BPDN_lee      (simInd, snrInd) = SE( SP.H, Hest_BPDN_lee,     SP, SP.SNR_lin_array(snrInd) );
        SP.SE_IHT_lee       (simInd, snrInd) = SE( SP.H, Hest_IHT_lee,      SP, SP.SNR_lin_array(snrInd) );
    end
    disp('')
end
end

