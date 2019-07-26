function [ A_t, A_r, beamPairs ] = ARM( Nt, Nr, Gt, Gr, angleType )
%ARM: Generate array response matrices.
switch angleType
    
    case 'spatial'
        %% Uniform angle grid in spatial angle (vartheta)
        tx_spatial_angle = (0:1:Gt-1)/Gt - 1/2;
        rx_spatial_angle = (0:1:Gr-1)/Gr - 1/2;
        tx_spatial_angle = [tx_spatial_angle(Gt/2+1:end) tx_spatial_angle(1:Gt/2)];
        rx_spatial_angle = [rx_spatial_angle(Gr/2+1:end) rx_spatial_angle(1:Gr/2)];
        % incident_angle_converted = acos(2*spatial_angle);
        A_t = zeros(Nt, Gt);
        A_r = zeros(Nr, Gr);
        for tx_angleInd = 1:length(tx_spatial_angle)
            at_temp = 1/sqrt(Nt) * exp(-1j*2*pi*(0:1:Nt-1)*tx_spatial_angle(tx_angleInd));
            A_t(:,tx_angleInd) = at_temp;
        end
        for rx_angleInd = 1:length(rx_spatial_angle)
            ar_temp = 1/sqrt(Nr) * exp(-1j*2*pi*(0:1:Nr-1)*rx_spatial_angle(rx_angleInd));
            A_r(:,rx_angleInd) = ar_temp;
        end
        
        
        tx_anglesInGrids = wrapTo2Pi(tx_spatial_angle*2*pi); % spatial angle (0 to 2pi)
        rx_anglesInGrids = wrapTo2Pi(rx_spatial_angle*2*pi);
        beamPairs = cell(Gt*Gr,1);
        ind = 0;
        for tx = 1:Gt
            for rx = 1:Gr
                ind = ind + 1; 
                beamPairs{ind} = [tx_anglesInGrids(tx), rx_anglesInGrids(rx)];
            end
        end
        disp('')
        
    case 'incident'
        %% Uniform angle grid in incident angle (theta)
        incident_angle = (0:1:G-1)/G * pi;
        A_t = zeros(Nt, G);
        A_r = zeros(Nr, G);
        for angleInd = 1:length(incident_angle);
            at_temp = 1/sqrt(Nt) * exp(-1j*pi*(0:1:Nt-1)*cos(incident_angle(angleInd)));
            ar_temp = 1/sqrt(Nr) * exp(-1j*pi*(0:1:Nr-1)*cos(incident_angle(angleInd)));
            A_t(:, angleInd) = at_temp;
            A_r(:, angleInd) = ar_temp;
        end
    otherwise
        error('Not proper angle type')
end

