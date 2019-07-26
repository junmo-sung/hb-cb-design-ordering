 function [Qbar, QbarNorm, Sigma, Q, A_t, A_r,  Wbar_diag, gridBeamPairs] = mtc_cb(SP)
[A_t, A_r, gridBeamPairs] = ARM(SP.Nt, SP.Nr, SP.Gt, SP.Gr, 'spatial'); % Array reponse matrices generation
NtBeam = SP.Mt;
NrBeam = SP.Mr*SP.Lr;
NtBlck = SP.Nt/SP.Lt;
NrBlck = SP.Nr/SP.Lr;
% X = eye(NtBeam);
switch SP.rfArchitecture
    case {'PS', 'lens'}
        Frf = 1/sqrt(SP.Nt) * dftmtx(SP.Nt);
        Wrf = 1/sqrt(SP.Nr) * dftmtx(SP.Nr);
        Frf = Frf(:, randperm(SP.Nt));
        Wrf = Wrf(:, randperm(SP.Nr));
        if SP.rfQuant == true && strcmp(SP.rfArchitecture, 'PS') % quantized angle
            Frf = quantizePSangles(SP, Frf);
            Wrf = quantizePSangles(SP, Wrf);
        end
        Fbb = sqrt(1/SP.Lt)*dftmtx(SP.Lt)* ...
            [eye(NtBeam/NtBlck), zeros(NtBeam/NtBlck, SP.Lt-NtBeam/NtBlck)]'*...
            sqrt(1/(NtBeam/NtBlck))*dftmtx(NtBeam/NtBlck)';
        Wbb = sqrt(1/SP.Lr)*dftmtx(SP.Lr)*...
            [eye(NrBeam/NrBlck), zeros(NrBeam/NrBlck, SP.Lr-NrBeam/NrBlck)]'*...
            sqrt(1/(NrBeam/NrBlck))*dftmtx(NrBeam/NrBlck)';
    case 'switches'
        Frf = eye(SP.Nt);
        Wrf = eye(SP.Nr);
        Frf = Frf(:, randperm(SP.Nt));
        Wrf = Wrf(:, randperm(SP.Nr));
        Ilt = eye(SP.Lt);
        Fbb = Ilt(:, 1:SP.sym);
        Wbb = eye(SP.Lr);
end
FbbKron = kron(eye(NtBlck), Fbb);
WbbKron = kron(eye(NrBlck), Wbb);
F = Frf*FbbKron;
W = Wrf*WbbKron;
Q = kron( F.', W' );
Qbar = Q*kron(conj(A_t), A_r);
% getMutCoh(Q)
Wbar_cell = mat2cell(W, SP.Nr, SP.Lr*ones(1, SP.Nr/SP.Lr));
Wbar_diag = blkdiag(Wbar_cell{:}); % This block diagonal matrix is for noise generation.
SigmaVec = zeros(size(Qbar,2),1);
for index = 1:size(Qbar,2)
    %     Sigma(index, index) = 1/norm(Phibar(:,index));
    SigmaVec(index) = 1/norm(Qbar(:,index));
    %     PhibarNorm(:, index) = Phibar(:, index)*Sigma(index,index); % <=> PhibarNorm = Phibar*Sigma;
end
Sigma = sparse(1:size(Qbar,2), 1:size(Qbar,2), SigmaVec);
QbarNorm = Qbar*Sigma;
end

