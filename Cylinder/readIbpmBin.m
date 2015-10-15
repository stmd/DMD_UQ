function grid = readIbpmBin(filename)
% Function to read .bin restart file from ibpm

% Open file
fid = fopen(filename,'r');
% Read grid info
nx = fread(fid,1,'int');
ny = fread(fid,1,'int');
ngrid = fread(fid,1,'int');
dx = fread(fid,1,'double');
x0 = fread(fid,1,'double');
y0 = fread(fid,1,'double');
numPoints = fread(fid,1,'int');
% Preallocations
q = zeros(nx,ny);
omega = zeros(ngrid,nx,ny);
FX = zeros(numPoints,1);
FY = zeros(numPoints,1);
numFluxes = 2*nx*ny + nx + ny;
% Read flux
disp('Reading flux...');
for lev=1:ngrid
    for qind=1:numFluxes
        q(lev,qind) = fread(fid,1,'double');
    end
end
% Read scalar omega
disp('Reading omega...');
for lev=1:ngrid
    for i=2:nx
        for j=2:ny
            omega(lev,i,j) = fread(fid,1,'double');
        end
    end
end
% Read BoundaryVector f
disp('Reading BoundaryVector...');
for i=1:numPoints
    FX(i) = fread(fid,1,'double');
    FY(i) = fread(fid,1,'double');
end
% Read timestep and time
disp('Reading timestep and time...');
timestep = fread(fid,1,'int');
time = fread(fid,1,'double');
% Close file
fclose(fid);
% Return values in struct
grid.nx = nx;
grid.ny = ny;
grid.ngrid = ngrid;
grid.dx = dx;
grid.x0 = x0;
grid.y0 = y0;
grid.numPoints = numPoints;
grid.q = q;
grid.omega = omega;
grid.FX = FX;
grid.FY = FY;
grid.numFluxes = numFluxes;
grid.timestep = timestep;
grid.time = time;

end