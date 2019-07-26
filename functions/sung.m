function [ Phibar, PhibarNorm, Sigma, Phi, A_t, A_r,  W_diag, gridBeamPairs] = sung2( SP )
% SP.Mt = SP.Nt/SP.Lt * SP.sym;
% SP.Mr = SP.Nr/SP.Lr;
% SP.M  = SP.Mt*SP.Mr;

%% BB
% Fbb and symbols
% I           = eye(SP.Lt);
% % % Xbar_cell   = cell(1,SP.Lt);
% x_cell = cell(SP.sym,1);
% for fInd    = 1:SP.sym
%     x_tilde         = I(:,fInd);
%     [x_cell{fInd}, ~, ~]    = BBprecoderGen(SP, x_tilde); % Fbbx
% %     Xbar_cell{fInd} = Frf*kron(eye(SP.Nt/SP.Lt), Fbbx);
% end

% X_temp = cell2mat(x_cell');
% X_temp = dctmtx(SP.Lt)';
% X_temp = 1/sqrt(SP.Lt)*dftmtx(SP.Lt);
% X_temp = X_temp(:, randperm(SP.sym));
% x_cell = mat2cell(X_temp, SP.Lt, ones(SP.Lt,1));

Wbb = eye(SP.Lr);
% Wbb = BBcombinerGen(SP);
fileName = [SP.rfArchitecture, '_Nt', num2str(SP.Nt), '_Lt', num2str(SP.Lt), '_sym', num2str(SP.sym)];
load(['ordered_symbol_precoder_set/', fileName]); % load x_cell and Frf

% Wbb = sqrt(1/SP.Lr)*dftmtx(SP.Lr);

% RF precoder/combiner and Array Response Matrix generation
[A_t, A_r, gridBeamPairs] = ARM(SP.Nt, SP.Nr, SP.Gt, SP.Gr, 'spatial'); % Array reponse matrices generation
% A_t = A_t(:, randperm(SP.Gt));
% A_r = A_r(:, randperm(SP.Gr));
switch SP.rfArchitecture
    case 'PS'
        switch SP.rfBFtype
            case 'DFT'
%                 Frf = 1/sqrt(SP.Nt) * dftmtx(SP.Nt); % unquantized RF precoder PS angles
                Wrf = 1/sqrt(SP.Nr) * dftmtx(SP.Nr); % unquantized RF combiner PS angles
%                 [~, Frf] = orderFrf(x_cell, Frf, SP, 'first');
%                 Frf = Frf(:, randperm(SP.Nt));
            case 'DCT'
                Frf = dctmtx(SP.Nt); % unquantized RF precoder PS angles
                Wrf = dctmtx(SP.Nr); % unquantized RF precoder PS angles
                % need to populate for
            case 'Hadamard'
                % Hadamard consists of +-1. Thus 2bit resolution is enough for it.
                Frf = 1/sqrt(SP.Nt) * hadamard(SP.Nt); % unquantized RF precoder PS angles
                Wrf = 1/sqrt(SP.Nr) * hadamard(SP.Nr); % unquantized RF combiner PS angles
            case 'GD'
%                 Frf = 1/sqrt(SP.Nt) * exp(2j*pi*rand(SP.Nt));
                  [Frf, x_cell] = loadXF(SP);
%                 Frf = getFRF_GD(SP, x_cell);
%                 [Frf, x_cell] = getFx_GD(SP, x_cell);
                Wrf = 1/sqrt(SP.Nr) * dftmtx(SP.Nr);
            otherwise
                error('Wrong analog beamformer type')
        end
        % quantized angle
        if SP.rfQuant == true
            Frf = quantizePSangles(SP, Frf);
            Wrf = quantizePSangles(SP, Wrf);
%             for m = 1:SP.Nt
%                 for n = 1:SP.Nt
%                     [~, Find] = min(abs(rfPossAngles - angle(Frf(m,n))));
%                     Frf(m,n) = 1/sqrt(SP.Nt)*exp(1j*rfPossAngles(Find));
%                 end
%             end
%             for m = 1:SP.Nr
%                 for n = 1:SP.Nr
%                     [~, Wind] = min(abs(rfPossAngles - angle(Wrf(m,n))));
%                     Wrf(m,n) = 1/sqrt(SP.Nr)*exp(1j*rfPossAngles(Wind));
%                 end
%             end
        end
        
    case 'switches'
%         Frf = eye(SP.Nt);
%         Frf = Frf(:, randperm(SP.Nt));
        Wrf = eye(SP.Nr);
    case 'lens'
%         Frf = 1/sqrt(SP.Nt)*dftmtx(SP.Nt);
        Wrf = 1/sqrt(SP.Nr)*dftmtx(SP.Nr);
    otherwise
        error('Wrong RF architecture')
end

%% BB
% x_mat = 1/sqrt(SP.Lt)*dftmtx(SP.Lt);
% x_cell = mat2cell(x_mat, SP.Lt, ones(SP.Lt,1));
% x_cell = getxCell_lee(SP);
% maxlatched = inf;
% for idx_iter = 1:SP.permIter
%     [x_cell_temp, maxval_temp] = getxCell(SP);
%     if maxval_temp < maxlatched
%         maxlatched = maxval_temp;
%         x_cell = x_cell_temp;
%     end
% end
% Wbb_rev = eye(SP.Lr);

%% Phi
% mutcoh = inf;

% for idx_iter = 1:SP.permIter
    % permute columns
%     testmat = zeros(SP.Nt);
%     Wrf = Wrf(:, randperm(SP.Nr));
    Frf_cell = mat2cell(Frf, SP.Nt, SP.Lt*ones(SP.Nt/SP.Lt, 1));
    Wrf_cell = mat2cell(Wrf, SP.Nr, SP.Lr*ones(SP.Nr/SP.Lr, 1));
%     x_mat = 1/sqrt(SP.Lt)*dftmtx(SP.Lt);
%     x_mat = eye(SP.Lt);
%     x_mat = x_mat(:, randperm(SP.Lt));
%     x_cell = mat2cell(x_mat, SP.Lt, ones(SP.Lt,1));
    Phi_cell = cell(SP.M, 1);
    W_cell = cell(SP.M, 1);
    s_cell = cell(1, SP.M);
    idx_Phi = 0;
    idx_x = 1;
    idx_Frf = 1;
    idx_Wrf = 0;
    
%     idx_x_remove = 0;
%     max_val = inf;
%     Unt = sqrt(1/SP.Nt) * dftmtx(SP.Nt);
%     Ult = sqrt(1/SP.Lt) * dftmtx(SP.Lt);
%     A = Unt * kron(eye(SP.Nt/SP.Lt), Ult);
%     Ilt = eye(SP.Lt);
%     for idx_Lt = 1:SP.Lt
%         B = kron( eye(SP.Nt/SP.Lt), diag( Ilt(:, idx_Lt) ) );
%         C = A*B*A';
%         max_val_current = max(max(abs(C)));
%         if max_val_current < max_val
%             max_val = max_val_current;
%             idx_x_remove = idx_Lt;
%         end
%     end
    
    while idx_Phi < SP.M
        idx_Phi = idx_Phi + 1;
        %     idx_x = idx_x + 1;
        %     if idx_x > SP.Lt, idx_x = 1; idx_Frf = idx_Frf+1; end
        %     if idx_Frf > SP.Nt/SP.Lt, idx_Frf = 1; idx_Wrf = idx_Wrf+1; end
        
        idx_Wrf = idx_Wrf + 1;
        if idx_Wrf > SP.Nr/SP.Lr, idx_Wrf = 1; idx_Frf = idx_Frf + 1; end
        if idx_Frf > SP.Nt/SP.Lt, idx_Frf = 1; idx_x = idx_x + 1; end
%         if idx_x == idx_x_remove, idx_x = idx_x + 1; end
%         idx_Frf = idx_Frf + 1;
%         if idx_Frf > SP.Nt/SP.Lt, idx_Frf = 1; idx_Wrf = idx_Wrf + 1; end
%         if idx_Wrf > SP.Nr/SP.Lr, idx_Wrf = 1; idx_x = idx_x + 1; end

        x   = x_cell{idx_x};
        Frf_m = Frf_cell{idx_Frf};
        Wrf_m = Wrf_cell{idx_Wrf};
        
        %     Frf = Frf_cell{ randi(SP.Nt/SP.Lt, 1) };
        %     x = x_cell{ randi(SP.Lt, 1) };
        %     Wrf = Wrf_cell{ randi(SP.Nr/SP.Lr, 1) };
        
        Phi_cell{idx_Phi} = kron((Frf_m*x).', (Wrf_m*Wbb)');
        W_cell{idx_Phi} = Wrf_m*Wbb;
        s_cell{idx_Phi} = Frf_m*x;
%         testmat = testmat + conj(Frf_m*x)*(Frf_m*x).';
    end
%     testmat = testmat / (SP.Nr/SP.Lr);
%     max(max((abs(testmat - diag(diag(testmat))))))
%     W_temp = zeros(SP.Nr);
%     for idx_m = 1:SP.M
%         W_temp = W_temp + (exp(-1j*5*pi/32))^(idx_m)*W_cell{idx_m}*W_cell{idx_m}';
%     end
%     
%     curr_mutcoh = getMutCoh( cell2mat(Phi_cell) )
%     mutcoh_vec(idx_iter) = curr_mutcoh;
%     if curr_mutcoh < mutcoh
%         mutcoh = curr_mutcoh;
%         Phi_cell_return = Phi_cell;
%         W_cell_return = W_cell;
%     end
% end
% Phi_cell = Phi_cell_return;
% W_cell = W_cell_return;

W_diag = blkdiag(W_cell{:});
Phi = cell2mat(Phi_cell);
% getMutCoh(Phi)
Phibar = Phi*kron(conj(A_t), A_r);

% Kt = Phibar*( diag(1./sqrt(diag(Phibar'*Phibar))) );
% KtKt = Kt'*Kt;
% max(max(abs(KtKt) - eye(SP.Gt*SP.Gr)))

% DELME4(SP, x_cell, Frf, A_t, A_r)
% [mu, nA] = getMutCoh(Phibar);

SigmaVec = zeros(size(Phibar,2),1);
for index = 1:size(Phibar,2)
    %     Sigma(index, index) = 1/norm(Phibar(:,index));
    SigmaVec(index) = 1/norm(Phibar(:,index));
    %     PhibarNorm(:, index) = Phibar(:, index)*Sigma(index,index); % <=> PhibarNorm = Phibar*Sigma;
end
Sigma = sparse(1:size(Phibar,2), 1:size(Phibar,2), SigmaVec);
PhibarNorm = Phibar*Sigma;

end

% %% BB precoder/combiner generation
% I           = eye(SP.Lt);
% Xbar_cell   = cell(1,SP.Lt);
% for fInd    = 1:SP.LtBar
%     x_tilde         = I(:,fInd);
%     [Fbbx, ~, ~]    = BBprecoderGen(SP, x_tilde);
%     Xbar_cell{fInd} = Frf*kron(eye(SP.Nt/SP.Lt), Fbbx);
% end
%
% Xbar    = cell2mat(Xbar_cell);
% Wbb     = BBcombinerGen(SP);
% W       = Wrf*kron(eye(SP.Nr/SP.Lr), Wbb);
% Wbar_cell = mat2cell(W, SP.Nr, SP.Lr*ones(1, SP.Nr/SP.Lr));
% Wbar_diag = blkdiag(Wbar_cell{:}); % This block diagonal matrix is for noise generation.
%
% %% Phi and Phi_bar Generation
% Phi    = kron(Xbar.', W'); % sensing matrix
% Phibar = Phi*kron(conj(A_t), A_r); % sensing*dictionary matrix (equivalent dictionary or eq. sens. mat.)
% PhibarNorm = zeros(size(Phibar));
% % Sigma  = zeros(size(Phibar,2));
% SigmaVec = zeros(size(Phibar,2),1);
% for index = 1:size(Phibar,2)
% %     Sigma(index, index) = 1/norm(Phibar(:,index));
%     SigmaVec(index) = 1/norm(Phibar(:,index));
% %     PhibarNorm(:, index) = Phibar(:, index)*Sigma(index,index); % <=> PhibarNorm = Phibar*Sigma;
% end
% Sigma = sparse(1:size(Phibar,2), 1:size(Phibar,2), SigmaVec);
% PhibarNorm = Phibar*Sigma;
% disp('')
%
%% Functions
function [Fbbx, Fbb, x] = BBprecoderGen(SP, x_input)
U = 1/sqrt(SP.Lt)*dftmtx(SP.Lt); V = U;
S_Fbb_vec   = ones(SP.Lt, 1).*x_input;
S_Fbb       = diag(S_Fbb_vec);
Fbbx        = U*S_Fbb*x_input;
Fbb         = U*S_Fbb*V';
x           = V*x_input;
end

function Wbb = BBcombinerGen(SP)
U = 1/sqrt(SP.Lr)*dftmtx(SP.Lr); V = U;
S_Wbb   = sqrt(SP.Nr/SP.Gr)*eye(SP.Lr);
Wbb     = U*S_Wbb*V';
end

