function matlab_pbsHandler(params,num)
% Matlab wrapper for job submission handling, to be called by Dakota

%------------------------------------------------------------------
% READ params.in (or params.in.num) from DAKOTA and set Matlab variables
%
% read params.in (no aprepro) -- just one param per line
% continuous design, then U.C. vars
% --> 1st cell of C has values, 2nd cell labels
% --> variables are in rows 2-->?
%------------------------------------------------------------------

paramsnum = [params '.' num2str(num)];
fid = fopen(paramsnum,'r');
C = textscan(fid,'%n%s');
fclose(fid);

% Set design variables -- could use vector notation
d = 4;
x = zeros(d,1);
for i=1:d
    x(i) = C{1}(i+1);
end
%------------------------------------------------------------------
% CALL your analysis code to get the function value
%------------------------------------------------------------------

newdir = pwd;
pbsfilename = [newdir '/pbsbase.dat'];
% Record input variables
UQparams = x';
save(strcat(newdir, '/params.out'), 'UQparams', '-ascii');
% Call function to compute DMD eigenvalue statistics
system(['echo "cd ' newdir  '" >> ' pbsfilename ]);
str = 'matlab -nodesktop -nodisplay -nosplash -r "calcDMDEvalStats(); exit"';
fid = fopen('pbsbase.dat','a');
fprintf(fid,'%s',str);
fprintf(fid,'\n');
% Write results to file
resultsfile = ['results.out.' num2str(num)];
fprintf(fid,'%s',['cp results.out ' resultsfile]);
fprintf(fid,'\n');
fclose(fid);
% Add to queue
system(['qsub ' pbsfilename]);


end
