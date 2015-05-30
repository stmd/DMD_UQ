#!/bin/csh -f


set num = `echo $argv[1] | cut -c 11-`
matlab -nojvm -nodesktop -nodisplay -nosplash -r "matlab_pbsHandler('params.in',$num); exit"
