function [ rateWF, rateEQ ] = SE( H, Hest, SP, SNR )

Ns = min(SP.Lt, SP.Lr);

if any(any( isnan(Hest) )) || any(any( isinf(Hest) ))
    rateEQ = -inf; % -inf indicates an error
    rateWF = -inf;
else
    [Uest, Sest, Vest] = svd(Hest);
    powerLevel = powerAllocation(Sest(:,1:Ns), Ns, SNR);
    Fest = Vest(:,1:Ns);
    West = Uest(:,1:Ns);
    Rn   = West'*West;
    Geq  = West'*H*Fest;
    Gwf  = West'*H*Fest*diag(powerLevel);
    
    rateEQ = abs(log2(det(eye(Ns) + SNR/Ns*inv(Rn)*(Geq*Geq'))));
    rateWF = abs(log2(det(eye(Ns) + SNR/Ns*inv(Rn)*(Gwf*Gwf'))));
    disp('')
    
end

    function powerLevel = powerAllocation(Sest, Ns, SNR)
        sigmaSqVec = diag(Sest).^2;
        powerLevel = -1;
        while min(powerLevel) < 0
            mu = Ns/length(sigmaSqVec) * (1 + 1/SNR * sum(1./sigmaSqVec));
            powerLevel = mu - Ns/SNR * 1./sigmaSqVec;
            sigmaSqVec = sigmaSqVec(1:end-1);
        end
        powerLevel = [powerLevel; zeros(Ns-length(powerLevel) ,1)];
        
        
    end

end

