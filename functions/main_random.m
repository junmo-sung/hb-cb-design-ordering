function [Hest_OMP, Hest_BPDN, Hest_IHT, Hest_LS, NMSE_LS, NMSE_OMP, NMSE_BPDN, NMSE_IHT] = main_random(SP, BF)

% [Phibar, Phi, A_t_bar, A_r_bar, Q_diag, SP] = md2016access(SP);
% SP.Mt = SP.Nt/SP.Lt * SP.Ns;
% SP.Mt = SP.Nt;
% SP.Mr = SP.Nr/SP.Lr;
% SP.M  = SP.Mt * SP.Mr;
% noisePow = SP.sigma^2 * SP.Mr*SP.Mt;
% noisePow = SP.Nr/SP.G*SP.sigma^2*SP.Lr*SP.Mr*SP.Mt;
noisePow = SP.sigma^2 * SP.M * SP.Lr;
AtAr = kron(conj(BF.AtBar), BF.ArBar);

switch SP.criteria
    case 'paths'
        criteria = SP.Np;
    case 'residual'
        criteria = sqrt(noisePow);
end

% for snrInd = 1:length(SP.SNR_db_array)
%     SP.rho = SP.SNR_lin_array(snrInd)*SP.sigma^2;
% Simulation iteration
%     for simInd = 1:SP.iter
%         fprintf('[mendez] SNR: %d (%d/%d), Iteration: %d/%d \n', SP.SNR_db_array(snrInd), snrInd, length(SP.SNR_db_array), simInd, SP.iter)
%         [H, A_t, A_r, AoD, AoA] = channelGen(SP);
%         noise = SP.sigma*sqrt(1/2)*...
%             (randn(SP.Mt*SP.Nr, 1) + 1j*randn(SP.Mt*SP.Nr, 1));
%         noise_matrix = SP.sigma*sqrt(1/2)*...
%             (randn(SP.Nr*SP.Mr, SP.Mt) + 1j*randn(SP.Nr*SP.Mr, SP.Mt));
% noise_combined = BF.Q_diag'*SP.noiseMatrix(:);
noise_vec_combined = BF.Q_diag'*SP.noiseVector;
h = SP.H(:);
y = sqrt(SP.rho)*BF.Phi*h + noise_vec_combined;

% LS
if ~isempty(find(strcmp(SP.algorithms, 'LS'), 1))
    tic;
    %         disp('Starting LS')
    h_hat_LS = 1/sqrt(SP.rho)*(BF.Phi\y);
    Hest_LS  = reshape(h_hat_LS, SP.Nr, SP.Nt);
    NMSE_LS = norm(h-h_hat_LS)^2/norm(h)^2;
    %         NMSE_LS(simInd, snrInd) = norm(h-h_hat_LS)^2/norm(h)^2;
    %         SP.NMSE_LS = NMSE_LS;
    %         10*log10(mean(NMSE_LS))
    elapsedTime = toc;
    fprintf('RND: LS took %.2f sec \n', elapsedTime);
else
    Hest_LS  = zeros(SP.Nr, SP.Nt);
    NMSE_LS  = 0;
end

% OMP
if ~isempty(find(strcmp(SP.algorithms, 'OMP'), 1))
    tic;
    %         disp('Starting OMP')
    [x_hat_OMP, ~, ~] = OMP(sqrt(SP.rho)*BF.Phibar, y, criteria); % OMP needs to be implemented by yourself.
    h_hat_OMP = AtAr*x_hat_OMP;
    Hest_OMP = reshape(h_hat_OMP, SP.Nr, SP.Nt);
    NMSE_OMP = norm(h-h_hat_OMP)^2/norm(h)^2;
    %         meanAngDiff = angleDiff(sparseToGridAngle(x_hat_OMP, BF.gridBeamPairs), SP.beamPairs);
    elapsedTime = toc;
    fprintf('RND: OMP took %.2f sec \n', elapsedTime);
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
    fprintf('RND: BPDN took %.2f sec \n', elapsedTime);
else
    Hest_BPDN = zeros(SP.Nr, SP.Nt);
    NMSE_BPDN = 0;
end

% IHT
if ~isempty(find(strcmp(SP.algorithms, 'IHT'), 1))
    tic;
    x_hat_IHT = DiagMatrixVectorProd(BF.Sigma, isht(BF.PhibarNorm, 20, y/sqrt(SP.rho), 1, criteria/sqrt(SP.rho), 10, 0));
    h_hat_IHT = AtAr*x_hat_IHT;
    Hest_IHT  = reshape(h_hat_IHT, SP.Nr, SP.Nt);
    NMSE_IHT  = norm(h-h_hat_IHT)^2/norm(h)^2;
    elapsedTime = toc;
    fprintf('RND: IHT took %.2f sec \n', elapsedTime);
else
    Hest_IHT = zeros(SP.Nr, SP.Nt);
    NMSE_IHT = 0;
end

% EM-GM-AMP
%         Hest_EMGMAMP  = zeros(SP.Nr, SP.Nt);
%         NMSE_EMGMAMP  = inf;
%         optEM.heavy_tailed = false;
%         x_hat_EMGMAMP   = EMGMAMP(y, sqrt(SP.rho)*BF.Phibar, optEM);
%         h_hat_EMGMAMP = AtAr*x_hat_EMGMAMP;
%         Hest_EMGMAMP  = reshape(h_hat_EMGMAMP, SP.Nr, SP.Nt);
%         NMSE_EMGMAMP = norm(h-h_hat_EMGMAMP)^2/norm(h)^2;



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

%         NMSE_OMP(simInd, snrInd) = norm(h-h_hat_OMP)^2/norm(h)^2;
%         SP.NMSE_OMP = NMSE_OMP;
%         10*log10(mean(NMSE_OMP))

%     end
% end
% SP.NMSE_LS  = NMSE_LS;
% SP.NMSE_OMP = NMSE_OMP;

