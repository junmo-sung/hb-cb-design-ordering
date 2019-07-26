function y = RC_shaping(T, beta, t)

y = zeros(size(t));

for i=1:length(t)
    if t(i) == T/(2*beta) || t(i) == -T/(2*beta)
        y(i) = pi/(4*T) * sinc(1/(2*beta));
    else
        y(i) = 1/T*sinc(t(i)/T)*cos(pi*beta*t(i)/T)/(1-(2*beta*t(i)/T)^2);
    end
end