%% LHS of the PCE surrogate

d = 4;
% Set parameter bounds
lower = [-1 -1 0 0]';
upper = [1 1 2 2]';
% Initialize LHS samples
samps = 1e5;
LHS = lhsdesign(samps,d);
X= zeros(samps,d);
for i=1:d
    X(:,i) = (upper(i)-lower(i))*LHS(:,i) + lower(i);
end
Y = evaluatePCESurrogate(X,lower,upper);

% Histogram of response statistics
%figure; hist(Y(:,1),50);
%figure; hist(Y(:,2),50);
%% Kernel-based smoothing of output PDFs

% Define kernel
Ker = @(x,mu,sig) exp(-0.5*((x-mu).^2)./(sig^2));
% Set output abscissa ranges
x1 = linspace(-0.2,1.2,1000)';
x2 = linspace(-0.15,0.25,1000)';
PDF1 = zeros(length(x1),1);
PDF2 = PDF1;
% Kernel smoothing
sig1 = 0.05;
sig2 = 0.01;
for i=1:length(Y)
    PDF1 = PDF1 + Ker(x1,Y(i,1),sig1);
    PDF2 = PDF2 + Ker(x2,Y(i,2),sig2);
    i
end
figure(10); plot(x1,PDF1./trapz(x1,PDF1),'b');
figure(11); plot(x2,PDF2./trapz(x2,PDF2),'b');

%% Calculate statistics of the surrogate

addpath('/home/adegenna/IcingShapeData/Dakota/templatedirREFINE/');
metric = CL;
% Pull out parts of data set of interest
p = [0:10:100]';
PCT = prctile(metric,p);
PCT = [p,PCT];
indLOW = find((metric>PCT(1,2)) & (metric<PCT(2,2)));
indHIGH = find((metric>PCT(10,2)) & (metric<PCT(11,2)));
low = randi(length(indLOW),20,1);
high = randi(length(indHIGH),20,1);
%LOW = X(indLOW(low),:);
%HIGH = X(indHIGH(high),:);
LOW = P(indLOW(low),:);
HIGH = P(indHIGH(high),:);
% Plot the corresponding airfoils
%
iceGoodX = zeros(350,20); iceGoodY = zeros(350,20);
iceBadX = zeros(350,20); iceBadY = zeros(350,20);
for i=1:20
    [airfoilInterp,~,~,~] = PODtoAirfoil(LOW(i,:));
    figure(19); hold on; plot(airfoilInterp(:,1),airfoilInterp(:,2),'Color','r');
    figure(19); xlim([-0.10,0.25]); axis equal;
    iceBadX(:,i) = airfoilInterp(:,1);
    iceBadY(:,i) = airfoilInterp(:,2);
end
for i=1:20
    [airfoilInterp,~,~,~] = PODtoAirfoil(HIGH(i,:));
    figure(20); hold on; plot(airfoilInterp(:,1),airfoilInterp(:,2),'Color','g');
    figure(20); xlim([-0.10,0.25]); axis equal;
    iceGoodX(:,i) = airfoilInterp(:,1);
    iceGoodY(:,i) = airfoilInterp(:,2);
    i
end
%}
%% Boxplot of POD coeffs

LOWall = P(indLOW,:);
HIGHall = P(indHIGH,:);
g1 = []; g2 = [];
for i=1:d
    LOWall(:,i) = 2*(LOWall(:,i)-lower(i))/(upper(i)-lower(i)) - 1;
    HIGHall(:,i) = 2*(HIGHall(:,i)-lower(i))/(upper(i)-lower(i)) - 1;
    g1 = [g1; i*ones(size(LOWall,1),1)];
    g2 = [g2; i*ones(size(HIGHall,1),1)];
end
figure; boxplot(LOWall(:),g1);
figure; boxplot(HIGHall(:),g2);

%% Estimate eigenvalue distribution in complex plane

d = 4;
% Set parameter bounds
lower = [-1 -1 0 0]';
upper = [1 1 2 2]';
% Initialize LHS samples of mean values
means = 1e4;
LHS = lhsdesign(means,d);
X= zeros(means,d);
for i=1:d
    X(:,i) = (upper(i)-lower(i))*LHS(:,i) + lower(i);
end
Y = evaluatePCESurrogate(X,lower,upper);
% Sample from each individual distribution
Ntrials = 100;
EIGS1 = [];
EIGS2 = [];
for i=1:means
    cov = [Y(i,3),Y(i,4);Y(i,4),Y(i,5)];
    [V,D] = eig(cov);
    eigs = diag(D);
    eigs(eigs<0) = 0;
    cov = V*diag(eigs)*V';
    samples = mvnrnd([Y(i,1);Y(i,2)],cov,Ntrials);
    EIGS1 = [EIGS1; samples];
    
    cov = [Y(i,8),Y(i,9);Y(i,9),Y(i,10)];
    [V,D] = eig(cov);
    eigs = diag(D);
    eigs(eigs<0) = 0;
    cov = V*diag(eigs)*V';
    samples = mvnrnd([Y(i,6);Y(i,7)],cov,Ntrials);
    EIGS2 = [EIGS2; samples];
    i
end

%% Compute actual eigenvalue locations at sparse grid locations

basedir = '/home/adegenna/SharedProj/DMD_UQ/DakotaFiles/workdir.';
Eigs1 = [];
for i=1:337
    directory = [basedir,num2str(i)];
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
    Eigs1 = [Eigs1; EigsTLS(1,:)'];
    i
end




