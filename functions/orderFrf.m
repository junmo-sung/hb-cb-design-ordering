function [x_chosen_cell_return, F] = orderFrf(x_cell, F, SP, opt)
% X = getMatSum(Xbar_cell);
% Xkron = kron(eye(SP.NtBlck), X);
Nt = SP.Nt;
I = eye(Nt);
% F_or(:,1) = F(:,1);
% F_un(:,1) = [];
if nargin < 4
    if strcmp(SP.rfArchitecture, 'PS')
        V_temp = nchoosek(2:SP.Lt, SP.sym-1);
        V = [ones(size(V_temp,1),1), V_temp];
    else
        V = nchoosek(1:SP.Lt, SP.sym);
    end
elseif strcmp(opt, 'first')
    V = 1:SP.sym;
else
    error('aa')
end

coherence_min1 = inf;
for idx_i = 1:size(V,1)
    fprintf('%d / %d \n', idx_i, size(V,1));
    x_chosen_cell = x_cell(V(idx_i,:));
    X = getMatSum(x_chosen_cell);
    Xkron = kron(eye(SP.NtBlck), X);
    
    % Column ordering
    F_un = F;
    F_or = zeros(size(F));
    for idx_n = 1:Nt
        F_test = F_or;
        coherence_min2 = inf;
        col_idx = 0;
        for idx_k = 1:size(F_un, 2)
            F_test(:,idx_n) = F_un(:,idx_k);
            S = F_test*Xkron*F_test';
            Sdia = diag(S);
            Sbis = diag(sqrt(1./Sdia));
            SS = Sbis*S*Sbis;
            coherence = max(max(abs(SS - I)));
            if coherence < coherence_min2
                coherence_min2 = coherence;
                col_idx = idx_k;
            end
        end
        if ~col_idx, break, end
        F_or(:,idx_n) = F_un(:,col_idx);
        F_un(:,col_idx) = [];
    end
    
    if coherence_min2 < coherence_min1
        coherence_min1 = coherence_min2;
        x_chosen_cell_return = x_cell(V(idx_i,:));
%         x_cell_idx = V(idx_i, :);
        F = F_or;
    end
end
% F = F_or;
% x_chosen_cell = x_cell(x_cell_idx);
end

function matsum = getMatSum(x_cell)
matsum = zeros(size(x_cell{1}, 1));
k = length(x_cell);
for idx_k = 1:k
    matsum = matsum + x_cell{idx_k}*x_cell{idx_k}';
end
end