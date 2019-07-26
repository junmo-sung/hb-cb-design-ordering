function [ C ] = khatriRao( A, B )
if size(A, 2) ~= size(B, 2)
    error('Input matrices have diffrent number of columns.')
end

M = size(A,1) * size(B,1);
N = size(A,2);
C = zeros(M,N);

for ind = 1:N
    C(:,ind) = kron(A(:,ind), B(:,ind));
end

end

