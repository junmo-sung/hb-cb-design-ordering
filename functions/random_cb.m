function [ Phibar, PhibarNorm, Sigma, Phi, A_t, A_r, W_diag, gridBeamPairs] = random_cb( SP )
rfPossAngles = (1:2^SP.rfQuantBits)*2*pi/2^SP.rfQuantBits - pi;
[A_t, A_r, gridBeamPairs] = ARM(SP.Nt, SP.Nr, SP.Gt, SP.Gr, 'spatial'); % Array reponse matrices generation

[Phi_cell, W_diag] = genPhiCell(SP, rfPossAngles);
Phi = cell2mat(Phi_cell);
Phibar = Phi*kron(conj(A_t), A_r);
SigmaVec = zeros(size(Phibar,2),1);
for index = 1:size(Phibar,2)
%     Sigma(index, index) = 1/norm(Phibar(:,index));
    SigmaVec(index) = 1/norm(Phibar(:,index));
%     PhibarNorm(:, index) = Phibar(:, index)*Sigma(index,index); % <=> PhibarNorm = Phibar*Sigma;
end
Sigma = sparse(1:size(Phibar,2), 1:size(Phibar,2), SigmaVec);
PhibarNorm = Phibar*Sigma;

end


function [Phi_cell, W_diag] = genPhiCell(SP, rfPossAngles)
Phi_cell = cell(SP.M, 1);
W_cell = cell(SP.M, 1);
switch SP.rfArchitecture
    case 'PS'
        for idx_m = 1:SP.M
            F = 1/sqrt(SP.Nt) * exp(1j*2*pi*rand(SP.Nt, SP.Lt));
            W = 1/sqrt(SP.Nr) * exp(1j*2*pi*rand(SP.Nr, SP.Lr));
            s = 1/sqrt(SP.Lt) * exp(1j*2*pi*rand(SP.Lt, 1));
            if SP.rfQuant
                F = quantPSangle(F, SP.Nt, rfPossAngles);
                W = quantPSangle(W, SP.Nr, rfPossAngles);
            end
            Phi_cell{idx_m} = kron( (F*s).', W' );
            W_cell{idx_m} = W;
        end
    case 'switches'
        Int = eye(SP.Nt);
        Inr = eye(SP.Nr);
        for idx_m = 1:SP.M
            Frand = Int(:, randperm(SP.Nt));
            Wrand = Inr(:, randperm(SP.Nr));
            F = Frand(:, 1:SP.Lt);
            W = Wrand(:, 1:SP.Lr);
%             s = 1/sqrt(SP.Lt) * exp(1j*2*pi*rand(SP.Lr,1));
            s = [1; zeros(SP.Lt-1,1)];
            Phi_cell{idx_m} = kron( (F*s).', W');
            W_cell{idx_m} = W;
        end
    case 'lens'
        Fdft = sqrt(1/SP.Nt)*dftmtx(SP.Nt);
        Wdft = sqrt(1/SP.Nr)*dftmtx(SP.Nr);
        for idx_m = 1:SP.M
            Frand = Fdft(:, randperm(SP.Nt));
            Wrand = Wdft(:, randperm(SP.Nr));
            F = Frand(:, 1:SP.Lt);
            W = Wrand(:, 1:SP.Lr);
            s = sqrt(1/SP.Lt) * exp(2j*pi*rand(SP.Lt,1));
            Phi_cell{idx_m} = kron( (F*s).', W');
            W_cell{idx_m} = W;
        end
        
end
W_diag = blkdiag(W_cell{:});
end

function M_quant = quantPSangle(M, mag, rfPossAngles)
[m,n] = size(M);
M_quant = zeros(size(M));
for idx_m = 1:m
    for idx_n = 1:n
        [~, Find] = min(abs(rfPossAngles - angle(M(idx_m,idx_n))));
        M_quant(idx_m,idx_n) = 1/sqrt(mag) * exp(1j*rfPossAngles(Find));
    end
end
end
