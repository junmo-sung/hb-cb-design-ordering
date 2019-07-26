close all
clear
clc

addpath functions

load_file = false;
SP.rfArchitecture   = 'PS'; 
SP.Nt = 64;
SP.Lt = 8;
SP.NtBlck = SP.Nt/SP.Lt;
SP.sym = 2;
SP.Mt = SP.NtBlck*SP.sym;
iteration = 3000;
fileName = [SP.rfArchitecture, '_Nt', num2str(SP.Nt), '_Lt', num2str(SP.Lt), '_sym', num2str(SP.sym)];

if ~load_file
    %% Ordering
    switch SP.rfArchitecture
        case 'PS'
            Frf = 1/sqrt(SP.Nt) * dftmtx(SP.Nt); % unquantized RF precoder PS angles
            X_temp = 1/sqrt(SP.Lt)*dftmtx(SP.Lt);
            x_cell = mat2cell(X_temp, SP.Lt, ones(SP.Lt,1));
            [x_cell, Frf] = orderFrf(x_cell, Frf, SP);
        case 'switches'
            Frf = eye(SP.Nt);
%             X_temp = 1/sqrt(SP.Lt)*dftmtx(SP.Lt);
%             X_temp = hadamard(SP.Lt);
            X_temp = eye(SP.Lt);
            x_cell = mat2cell(X_temp, SP.Lt, ones(SP.Lt,1));
            [x_cell, Frf] = orderFrf(x_cell, Frf, SP);
            
%             X_temp = eye(SP.Lt);
%             X_temp = X_temp(:, 1:1:4);
%             x_cell = mat2cell(X_temp, SP.Lt, ones(SP.sym,1));
%             X_temp = ones(SP.Lt);
%             

    end
    
    save(['temp_ordered_symbol_precoder_set/', fileName], 'x_cell', 'Frf');
else
    load(['temp_ordered_symbol_precoder_set/', fileName]);
end
X = getMatSum(x_cell);
coherence_proposed = getSimpleCoherence(SP, X, Frf)


%% Random permutation
% X_temp = 1/sqrt(SP.Lt)*dftmtx(SP.Lt);
% x_cell = mat2cell(X_temp, SP.Lt, ones(SP.Lt,1));
switch SP.rfArchitecture
    case 'PS'
        Frf = 1/sqrt(SP.Nt) * dftmtx(SP.Nt);
end
Fbb = sqrt(1/SP.Lt)*dftmtx(SP.Lt)* ...
            [eye(SP.Mt/SP.NtBlck), zeros(SP.Mt/SP.NtBlck, SP.Lt-SP.Mt/SP.NtBlck)]'*...
            sqrt(1/(SP.Mt/SP.NtBlck))*dftmtx(SP.Mt/SP.NtBlck)';
x_cell = mat2cell(Fbb, SP.Lt, ones(SP.sym,1));
X = getMatSum(x_cell);
% V = nchoosek(1:SP.Lt, SP.sym);
coherence_vec = zeros(iteration,1);
for idx = 1:iteration
    disp(idx)
%     x_chosen_cell = x_cell(
    Frf_perm = Frf(:, randperm(SP.Nt));
    coherence_vec(idx) = getSimpleCoherence(SP, X, Frf_perm);
end

%% Plotting
close all
figureSize = [8,5]; %inches
figure
hold on
set(gca, 'FontSize', 13)
set(gcf, 'Units', 'inches')
pos = get(gcf, 'position');
pos(3:4) = figureSize;
set(gcf, 'position', pos)
[counts, centers] = hist(coherence_vec, 100);
prob = counts/iteration;

bar(centers, prob)

[~, closest_index] = min(abs(centers - mean(coherence_vec)));
nearest_coherence = centers(closest_index);
nearest_coherence_prob = prob(closest_index);
plot(nearest_coherence, nearest_coherence_prob, 'ks', 'MarkerSize', 13, 'MarkerFaceColor', 'g')

[~, closest_index] = min(abs(centers - coherence_proposed));
nearest_coherence = centers(closest_index);
nearest_coherence_prob = prob(closest_index);
plot(nearest_coherence, nearest_coherence_prob, 'ko', 'MarkerSize', 13, 'MarkerFaceColor', 'r')

hold off
grid on
xlabel('Mutual Coherence')
ylabel('Probability')
% title('Nt64 / Nr16 / Lt8 / Lr8 / G64 / Np4 / Q16')
legend('Distribution of random permutation', 'Mean of random permutation', 'Proposed algorithm')

coherence_proposed
mean(coherence_vec)
%% Functions

function coherence = getSimpleCoherence(SP, X, F)
Xkron = kron(eye(SP.NtBlck), X);
S = conj(F)*Xkron*F.';
Sdia = diag(S);
Sbis = diag(sqrt(1./Sdia));
SS = Sbis*S*Sbis;
coherence = max(max(abs(SS - eye(SP.Nt))));
end

function matsum = getMatSum(x_cell)
matsum = zeros(size(x_cell{1}, 1));
k = length(x_cell);
for idx_k = 1:k
    matsum = matsum + x_cell{idx_k}*x_cell{idx_k}';
end
end