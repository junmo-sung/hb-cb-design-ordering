function M = quantizePSangles(SP, M)
rfPossAngles = (1:2^SP.rfQuantBits)*2*pi/2^SP.rfQuantBits - pi;
N = size(M,1);
for m = 1:N
    for n = 1:N
        [~, Find] = min(abs(rfPossAngles - angle(M(m,n))));
        M(m,n) = 1/sqrt(N)*exp(1j*rfPossAngles(Find));
    end
end
% for m = 1:SP.Nr
%     for n = 1:SP.Nr
%         [~, Wind] = min(abs(rfPossAngles - angle(Wrf(m,n))));
%         Wrf(m,n) = 1/sqrt(SP.Nr)*exp(1j*rfPossAngles(Wind));
%     end
% end
end

