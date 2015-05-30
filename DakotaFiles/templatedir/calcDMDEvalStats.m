function calcDMDEvalStats()
% (1) Perform Monte Carlo trials on DMD data with noise (specified by
% 'directory') added
% (2) Compute empirical mean and covariance of DMD eigenvalue distributions
% (3) Write means/covariances to file

directory = pwd;
% Set up dynamical system
Ac = [1 -2; 1 -1];
dt = 0.1;
A = expm(Ac*dt);
x0 = [ 1; 0.1];
Nsteps = 100;
x = zeros(2,Nsteps);
x(:,1) = x0;
% Attain noise-free data
for kk = 1:(Nsteps-1)
    x(:,kk+1) = A*x(:,kk);
end
% Get bias and scaling for noise
infile = strcat(directory,'/params.out');
X = importdata(infile);
% Monte Carlo trials
r = 2;
Ntrials = 500;
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
    % Separate eigenvalues into two groups
    eval = eig(Admdtls);
    [~,i1] = min(abs(eval-1i));
    [~,i2] = min(abs(eval+1i));
    EigsTLS(:,nn) = [eval(i1); eval(i2)];
end
% Compute statistics of eigenvalue distributions
real1 = real(EigsTLS(1,:));
imag1 = imag(EigsTLS(1,:));
real2 = real(EigsTLS(2,:));
imag2 = imag(EigsTLS(2,:));
MU1(1) = mean(real1);
MU1(2) = mean(imag1);
MU2(1) = mean(real2);
MU2(2) = mean(imag2);
COV1 = cov(real1,imag1);
CXX(1) = COV1(1,1);
CXY(1) = COV1(1,2);
CYY(1) = COV1(2,2);
COV2 = cov(real2,imag2);
CXX(2) = COV2(1,1);
CXY(2) = COV2(1,2);
CYY(2) = COV2(2,2);
% Write to file
outfile = strcat(directory,'/results.out');
fid=fopen(outfile,'w');
fprintf(fid,'%f MeanReal1 %f MeanImag1 %f Cxx1 %f Cxy1 %f Cyy1 %f MeanReal2 %f MeanImag2 %f Cxx2 %f Cxy2 %f Cyy2',...
    MU1(1),MU1(2),CXX(1),CXY(1),CYY(1),MU2(1),MU2(2),CXX(2),CXY(2),CYY(2));
fclose(fid);
% Plot
figure(1); plot(real1,imag1,'b.',real2,imag2,'r.');



end