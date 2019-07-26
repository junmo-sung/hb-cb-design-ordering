function [Hest_OMP, Hest_BPDN, Hest_IHT, Hest_LS, H_hat_or, NMSE_LS, NMSE_OMP, NMSE_BPDN, NMSE_IHT, NMSE_OR] = main_mtc(SP, BF)

% [Phibar, Phi, A_t_bar, A_r_bar, Wbar_diag] = sung2(SP);
noisePow = SP.sigma^2*SP.Lr*SP.Mr*SP.Mt;
% noisePow = SP.sigma^2 * SP.M * SP.Lr;
AtAr = kron(conj(BF.AtBar), BF.ArBar);

switch SP.criteria
    case 'paths'
        criteria = SP.Np;
    case 'residual'
        criteria = sqrt(noisePow);
end

% SP.Mt = SP.Nt/SP.Lt * SP.Ns;
% SP.Mt = SP.Nt;
% SP.Mr = SP.Nr/SP.Lr;
% SP.M  = SP.Mt * SP.Mr;

% switch SP.criteria
%     case 'paths'
%         criteria = SP.Np;
%     case 'residual'
%         criteria = sqrt(noisePow);
% end

% for snrInd = 1:length(SP.SNR_db_array)
%     SP.rho = SP.SNR_lin_array(snrInd)*SP.sigma^2;
%     % Simulation iteration
%     for simInd = 1:SP.iter
%         fprintf('[sung] SNR: %d (%d/%d), Iteration: %d/%d \n', SP.SNR_db_array(snrInd), snrInd, length(SP.SNR_db_array), simInd, SP.iter)
%         y_mat = zeros(SP.L_r*SP.N, SP.M);
%         y_mat = mat2cell(y_mat, SP.L_r*SP.N, ones(1,SP.M));
%         Phi = zeros(SP.N*SP.M*SP.L_r, SP.N_c*SP.N_r*SP.N_t);
%         Phi = mat2cell(Phi, SP.N*SP.L_r*ones(1, SP.M), SP.N_c*SP.N_r*SP.N_t);
%         receivedSigLen = SP.L_r*SP.N;
%         PhiSize = SP.N*SP.L_r;
%         [H, A_t, A_r, AoD, AoA] = channelGen(SP);
% Frame iteration
%         for frameInd = 1:SP.M
%% Receive vector generation
%         [F, precoderPhase] = precoderGen(SP);
%         [W, combinerPhase] = combinerGen(SP);
%         S = trainSigMatGen(SP.N, SP.N_c, SP.N_s, SP.trainType);
%         [y, y_withoutNoise] = receivedSigVec(SP, H, F, W, S);
%         y_mat{frameInd} = y;

%         noise_matrix = SP.sigma*sqrt(1/2)*...
%             (randn(SP.Nr*SP.Mr, SP.Mt) + 1j*randn(SP.Nr*SP.Mr, SP.Mt));
noiseMatrix = reshape(SP.noiseVector, SP.Nr*SP.Mr, SP.Mt);
noise_combined = BF.WBarDiag'*noiseMatrix;
noise_vec_combined = noise_combined(:);
% noise_vec_combined = BF.WBarDiag'*SP.noiseVector;

%         noise_vec_combined'*noise_vec_combined
%         norm(noise_combined)^2

h = SP.H(:);
y = sqrt(SP.rho)*BF.Phi*h + noise_vec_combined;
%% Measurement matrix generation
%         Phi_frame = kron(S*kron(eye(SP.N_c), F.'), W');
%         Phi{frameInd} = ...
%             sqrt(SP.rho)*Phi_frame;
%         %         end
%         y_mat = cell2mat(y_mat);
%         y_vec = y_mat(:);
%         y_quant = quantization(y_vec, SP.ADCbits, max(abs([real(y_vec);imag(y_vec)])));
%         Phi = cell2mat(Phi);
%         h = H(:);
%         PhiPsi = Phi*Psi;

% LS
if ~isempty(find(strcmp(SP.algorithms, 'LS'), 1))
    tic;
    %         disp('Starting LS')
    h_hat_LS = 1/sqrt(SP.rho)*(BF.Phi\y);
    Hest_LS  = reshape(h_hat_LS, SP.Nr, SP.Nt);
    %         NMSE_LS(simInd, snrInd) = norm(h-h_hat_LS)^2/norm(h)^2;
    NMSE_LS = norm(h-h_hat_LS)^2/norm(h)^2;
    %         SP.NMSE_LS = NMSE_LS;
    %         10*log10(mean(NMSE_LS))
    elapsedTime = toc;
    fprintf('DET: LS took %.2f sec \n', elapsedTime);
else
    Hest_LS  = zeros(SP.Nr, SP.Nt);
    NMSE_LS  = 0;
end

% OMP
if ~isempty(find(strcmp(SP.algorithms, 'OMP'), 1))
    %         disp('Starting OMP')
    tic;
%     opts.printEvery = 1e5;
%     opts.verbose    = 0;
    [x_hat_OMP, ~, ~] = OMP(sqrt(SP.rho)*BF.Phibar, y, criteria); % OMP needs to be implemented by yourself.
    h_hat_OMP = AtAr*x_hat_OMP;
    Hest_OMP = reshape(h_hat_OMP, SP.Nr, SP.Nt);
    %         NMSE_OMP(simInd, snrInd) = norm(h-h_hat_OMP)^2/norm(h)^2;
    NMSE_OMP = norm(h-h_hat_OMP)^2/norm(h)^2;
    % meanAngDiff = angleDiff(sparseToGridAngle(x_hat_OMP, BF.gridBeamPairs), SP.beamPairs);
    elapsedTime = toc;
    fprintf('DET: OMP took %.2f sec \n', elapsedTime);
else
    Hest_OMP = zeros(SP.Nr, SP.Nt);
    NMSE_OMP = 0;
end

% BPDN (http://www.cs.ubc.ca/~mpf/spgl1/download.html)
if ~isempty(find(strcmp(SP.algorithms, 'BPDN'), 1))
    tic;
    opts = spgSetParms('verbosity',0);
    x_hat_BPDN = DiagMatrixVectorProd(BF.Sigma, spg_bpdn(BF.PhibarNorm, y/sqrt(SP.rho), criteria/sqrt(SP.rho), opts));
    h_hat_BPDN = AtAr*x_hat_BPDN;
    Hest_BPDN  = reshape(h_hat_BPDN, SP.Nr, SP.Nt);
    NMSE_BPDN  = norm(h-h_hat_BPDN)^2/norm(h)^2;
    elapsedTime = toc;
    fprintf('DET: BPDN took %.2f sec \n', elapsedTime);
else
    Hest_BPDN = zeros(SP.Nr, SP.Nt);
    NMSE_BPDN = 0;
    
end
% IHT
if ~isempty(find(strcmp(SP.algorithms, 'IHT'), 1))
    tic;
    % x_hat_init = randn(size(BF.Phibar,2), 1);
    % x_hat_IHT = IHT(@(x)norm( sqrt(SP.rho)*BF.Phibar*x - y )^2, @(x)2*BF.Phibar'*(sqrt(SP.rho)*BF.Phibar-y), SP.Np, 2*max(eig(BF.Phibar'*BF.Phibar))+0.1, x_hat_init, 10);
    % options.T    =  criteria/sqrt(SP.rho);
    % options.Tmin = 4;
    % options.Tmax = 40;
    % options.thresh = 'hard';
    % x_hat_IHT = perform_iterative_thresholding(sqrt(SP.rho)*BF.Phibar, y, options);
    % x_hat_IHT = BF.Sigma*perform_iterative_thresholding(BF.PhibarNorm, y/sqrt(SP.rho), options);
    % sol_len = size(BF.Phibar, 2);
    % Psi = BF.PhibarNorm; lambda = 1; method = 1; xinitial = randn(sol_len, 1) + 1j*randn(sol_len,1); maxIt = 100; tol=criteria/sqrt(SP.rho); mun=size(BF.Phibar,1)/2;
    % x_hat_IHT = BF.Sigma*isht(Psi, Psi', y/sqrt(SP.rho), lambda, method, xinitial, maxIt, tol, mun);
    x_hat_IHT = DiagMatrixVectorProd(BF.Sigma, isht(BF.PhibarNorm, 20, y/sqrt(SP.rho), 1, criteria/sqrt(SP.rho), 10, 0));
    h_hat_IHT = AtAr*x_hat_IHT;
    Hest_IHT  = reshape(h_hat_IHT, SP.Nr, SP.Nt);
    NMSE_IHT  = norm(h-h_hat_IHT)^2/norm(h)^2;
    % close all; figure; plot(abs(x_hat_OMP)); figure; plot(abs(x_hat_IHT));
    elapsedTime = toc;
    fprintf('DET: IHT took %.2f sec \n', elapsedTime);
    disp('')
else
    Hest_IHT = zeros(SP.Nr, SP.Nt);
    NMSE_IHT = 0;
end


% GAMP
% Hest_EMGMAMP  = zeros(SP.Nr, SP.Nt);
% NMSE_EMGMAMP  = inf;
% % noise_variance = var(noise_vec_combined);
% % noise_variance = var(SP.noiseMatrix(:));
% %     x_hat_GAMP = CGAMP(sqrt(SP.rho)*BF.Phibar, y, 0, mean(diag(h*h')), SP.Np/SP.G^2, noise_variance);
% optEM.heavy_tailed = false;
% x_hat_EMGMAMP   = EMGMAMP(y, sqrt(SP.rho)*BF.Phibar, optEM);
% % h_hat_GAMP = AtAr*x_hat_GAMP;
% h_hat_EMGMAMP = AtAr*x_hat_EMGMAMP;
% % Hest_GAMP  = reshape(h_hat_GAMP, SP.Nr, SP.Nt);
% % NMSE_GAMP  = norm(h-h_hat_GAMP)^2/norm(h)^2;
% Hest_EMGMAMP = reshape(h_hat_EMGMAMP, SP.Nr, SP.Nt);
% NMSE_EMGMAMP = norm(h-h_hat_EMGMAMP)^2/norm(h)^2;

% Oracle
PhiOr = BF.Phi*khatriRao(conj(SP.A_t), SP.A_r);
h_hat_or = 1/sqrt(SP.rho)*(PhiOr\y);
H_hat_or = SP.A_r*diag(h_hat_or)*SP.A_t';
NMSE_OR = norm(h-H_hat_or(:))^2 / norm(h)^2;

end

function meanDiff = angleDiff(angleCellA, angleCellB) % A: estimation, B: true
%     diffCell = cell(length(angleCellA),1);
    diffVec = zeros(length(angleCellA), 1);
    for indA = 1:length(angleCellA)
        minDiff = [inf,inf];
        for indB = 1:length(angleCellB)
            diff = angleCellB{indB} - angleCellA{indA};
            if norm(diff) < norm(minDiff)
                minDiff = diff;
            end
        end
%         diffCell{indA} = minDiff;
        diffVec(indA) = (norm(minDiff))^2;
    end
    meanDiff = mean(diffVec);
end

function angleCell = sparseToGridAngle(sparseVec, gridAngleCell)
    sparseIndVec = find(sparseVec ~= 0);
    angleCell = cell(length(sparseIndVec), 1);
    for ind = 1:length(sparseIndVec)
        angleCell{ind} = gridAngleCell{sparseIndVec(ind)};
    end
end

function product = DiagMatrixVectorProd(diagonal_matrix, vector)
    product = zeros(size(vector));
    for ind = 1:length(vector)
        product(ind) = diagonal_matrix(ind,ind) * vector(ind);
    end
end
%         SP.NMSE_OMP = NMSE_OMP;
%         10*log10(mean(NMSE_OMP))
%     end
% end
% SP.NMSE_LS  = NMSE_LS;
% SP.NMSE_OMP = NMSE_OMP;