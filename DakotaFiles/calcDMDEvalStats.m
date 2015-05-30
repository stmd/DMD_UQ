function calcDMDEvalStats(directory)
% (1) Perform Monte Carlo trials on DMD data with noise (specified by
% 'directory') added
% (2) Compute empirical mean and covariance of DMD eigenvalue distributions
% (3) Write means/covariances to file

% Set up dynamical system
Ac = [1 -2; 1 -1];
dt = 0.1;
A = expm(Ac*dt);
EigsTrue = eig(A);
x0 = [ 1; 0.1];
Nsteps = 100;
x = zeros(2,Nsteps);
x(:,1) = x0;
% Attain noise-free data
for kk = 1:(Nsteps-1)
    x(:,kk+1) = A*x(:,kk);
end
% Get bias and scaling for noise
X = importdata(directory);
% Monte Carlo trials
r = 2;
Ntrials = 100;
s = 0.2; %noise level
EigsTLS = zeros(2,Ntrials);
for nn = 1:Ntrials
    Noise = s*randn(size(x));
    % Bias/scale the noise
    Noise(1,:) = X(1) + X(3)*Noise(1,:);
    Noise(2,:) = X(2) + X(4)*Noise(2,:);
    xn = x + Noise;
    % Compute DMD eigenvalues
    [Phin, Lambdan, U, S, V, Atilden] = DMDext(xn,r);
    Admdtls  = dmd_tls(xn);
    EigsTLS(:,nn) = sort(eig(Admdtls));
end



end