function Y = evaluatePCESurrogate(X,lower,upper)
% Function to evaluate PCE surrogate

directory = '/home/adegenna/SharedProj/DMD_UQ/DakotaFiles/TOL1em4';
outputs = 10;
% Import PCE coefficients
DATA = importdata([directory,'/PCEOutputCoeffs.dat'],' ');
COEFFS = DATA(:,1:outputs);
IND = DATA(:,outputs+1:end);
% Evaluate surrogate at sample locations
Y = zeros(length(X),outputs);
for i=1:length(COEFFS)
    phi = ones(length(X),1);
    for j=1:length(lower)
        %norm = 2/(2*IND(j)+1);
        XI = (2/(upper(j)-lower(j)))*X(:,j) - (upper(j)+lower(j))/(upper(j)-lower(j));
        phi = phi.*LEGENDRE(IND(i,j),XI)';
    end
    for j=1:outputs
        Y(:,j) = Y(:,j) + COEFFS(i,j)*phi;
    end
    i
end


end

