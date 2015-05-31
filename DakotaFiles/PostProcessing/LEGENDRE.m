function evals = LEGENDRE(N,XI)
% N denotes the Nth Legendre polynomial (counting starting at zero; highest
% order monomial term is order N)
% XI is the point for evaluation of the polynomial

K = zeros(N+1,1);
% Expansion is L_N(XI) = K(1)*XI^N + K(2)*XI^(N-1) + ... + K(N)*XI + K(N+1)
% Find expansion coefficients (NOT VECTORIZED)
%{
count = 1; COUNT = [];
for i=0:floor(N/2)
    COUNT = [COUNT; i,count];
    KTMP(count,1) = (1/2^N)*(-1)^i*nchoosek(N,i)*nchoosek(2*(N-i),N);
    count = count + 1;
end
K(1:2:2*(count-1)-1) = KTMP;
%}

% Find expansion coefficients (VECTORIZED)
i = [0:floor(N/2)]';
count = [1:floor(N/2)+1]';
KTMP(:,1) = (1/2^N)*(-1).^i.*NCHOOSEK(N,i).*NCHOOSEK(2*(N-i),N);
K(1:2:2*(floor(N/2)+1)-1) = KTMP;

%{
count = 1;
for i=N:-1:0
    %MONOM(count) = XI^i;
    MONOM(count) = XI.^i;
    count = count + 1;
end
%}
i = [N:-1:0]';
[XXI,II] = meshgrid(XI,i);
MONOM = XXI.^II;

evals = K'*MONOM;

end