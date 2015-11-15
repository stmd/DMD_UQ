function Eigs_real_imag = compute_SL_Evals(params)
% Function to compute the eigenvalues of the Stuart-Landau equation
% INPUTS:
%   params: [mu,beta] (governing parameters)
% OUTPUTS:
%   Eigs_SL:   [real(Eig); imag(Eig)]
%              ordered as [real(lambda_1,...,lambda_n);
%              imag(lambda_1,...,lambda_n)]

% Extract parameter values
mu = params(1);
beta = params(2);
% Set time vector and other parameters
dt = 0.01;
t= 0:dt:10;
r0 = 2;
th0 = 0;
gamma = 1;
% Compute SL dynamics
r2 = r0*(mu./(r0^2+(mu-r0^2)*exp(-2*mu*t))).^0.5;
th = th0 + (gamma-mu*beta)*t+beta*log(r2/r0);
% Set up data matrix
minr = -1;
maxr = 0;
minn = -3;
maxn = 3;
Nrows = (maxn-minn+1)*(maxr-minr+1);
order = Nrows;
Y = zeros(Nrows,length(t));
k= 1;
for rr = minr:maxr
    for nn = minn:maxn
        Y(k,:) = r2.^(2*rr).*exp(1i*nn*th);
        k=k+1;
    end
end
% Compute eigenvalues
[~,LambdaE,~,~,~] = DMDext(Y,order);
Eigs = diag(LambdaE);
% Convert discrete --> continuous
Eigsc = log(Eigs)/dt;
% Sort eigenvalues
[~,EigsSortRealInd] = sort(real(Eigsc));
EigsSortReal = Eigs(EigsSortRealInd);
[~,EigsSortImag1Ind] = sort(imag(EigsSortReal(1:end/2)));
[~,EigsSortImag2Ind] = sort(imag(EigsSortReal((end/2+1):end)));
Eigs_SL = [EigsSortReal(EigsSortImag1Ind);EigsSortReal(EigsSortImag2Ind+(maxn-minn+1))];
% Split into real/imag parts
Eigs_real_imag(:,1) = real(Eigs_SL);
Eigs_real_imag(:,2) = imag(Eigs_SL);

end