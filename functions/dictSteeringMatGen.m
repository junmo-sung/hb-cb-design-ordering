function [A_tx, A_rx, DictAoD, DictAoA] = dictSteeringMatGen(N_t, N_r, G_t, G_r, virtualChannel)

if strcmp(virtualChannel, 'physical')
    a_t = @(N_t, aod) 1/sqrt(N_t)*exp(1j*pi*(0:N_t-1).'*cos(aod));
    a_r = @(N_r, aoa) 1/sqrt(N_r)*exp(1j*pi*(0:N_r-1).'*cos(aoa));
elseif strcmp(virtualChannel, 'spatial')
    a_t = @(N_t, aod) 1/sqrt(N_t)*exp(-1j*(0:N_t-1).'*aod);
    a_r = @(N_r, aoa) 1/sqrt(N_r)*exp(-1j*(0:N_r-1).'*aoa);
end

DictAoD = (0:G_t-1)*2*pi/G_t;
DictAoA = (0:G_r-1)*2*pi/G_r;

A_tx = a_t(N_t, DictAoD);
A_rx = a_r(N_r, DictAoA);