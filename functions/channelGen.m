function [H, A_t, A_r, beamPairs] = channelGen(SP)

if strcmp(SP.virtualChannel, 'physical')
    a_t = @(N_t, aod) 1/sqrt(N_t)*exp(1j*pi*(0:N_t-1).'*cos(aod));
    a_r = @(N_r, aoa) 1/sqrt(N_r)*exp(1j*pi*(0:N_r-1).'*cos(aoa));
elseif strcmp(SP.virtualChannel, 'spatial')
    a_t = @(N_t, aod) 1/sqrt(N_t)*exp(-1j*(0:N_t-1).'*aod);
    a_r = @(N_r, aoa) 1/sqrt(N_r)*exp(-1j*(0:N_r-1).'*aoa);
end
[~, ~, DictAoD, DictAoA] = dictSteeringMatGen(SP.Nt, SP.Nr, SP.Gt, SP.Gr, SP.virtualChannel);
if SP.integerTapDelay
    tau = randi(SP.Nc, SP.Np, 1)-1;
else
    tau = (SP.Nc-1)*rand(SP.Np, 1); % zeros(Np, 1);
end
alpha = sqrt(1/2)*(randn(SP.Np,1) + 1j*randn(SP.Np, 1)); % channel gain

if SP.BeamAlign
    AoDindex = randi(SP.Gt, SP.Np, 1);
    AoAindex = randi(SP.Gr, SP.Np, 1);
    AoD = DictAoD(AoDindex);
    AoA = DictAoA(AoAindex);
else
    AoD = 2*pi*rand(1, SP.Np);
    AoA = 2*pi*rand(1, SP.Np);
    beamPairs = {};
        
end

beamPairs = cell(SP.Np, 1);
for i = 1:SP.Np
    beamPairs{i} = [AoD(i), AoA(i)];
end

A_t = a_t(SP.Nt, AoD);
A_r = a_r(SP.Nr, AoA);
    
H = zeros(SP.Nr, SP.Nt, SP.Nc);
for d = 1:SP.Nc
    p     = RC_shaping(1, SP.beta, (d-1)-tau);
    Delta = diag(p.*alpha);
    H_temp = A_r*Delta*A_t';
    H(:,:,d) = sqrt(SP.Nr*SP.Nt/SP.Np)*H_temp; %/norm(H_temp, 'fro');
    disp('')
end
