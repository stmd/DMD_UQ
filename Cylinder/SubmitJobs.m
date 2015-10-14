% Submit batch jobs; read Reynolds number from file

Re = [
   97.6545
   88.4617
   61.5383
   52.3455];
% Define base/batch directories
basedir = '/home/adegenna/ibpmcontrol/ibpmcontrol';
batchdir = '/home/adegenna/ibpmcontrol/ibpmcontrol/output/Cylinder';
%
for i=1:length(Re)
    workdir = ['/home/adegenna/ibpmcontrol/ibpmcontrol/output/Cylinder/CylInitRe',num2str(Re(i))];
    cd(batchdir);
    mkdir(workdir);
    system(['cp pbsbase.dat ' workdir]);
    % qsub
    cd(basedir);
    system(['qsub -v OUTDIR=' workdir ',RE=' num2str(Re(i)) ' ' batchdir '/pbsbase.dat']);
end
%}