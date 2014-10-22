# 
#$ -cwd 
#$ -j y 
#$ -pe orte 12
#$ -S /bin/bash 
#
# Name the job #$ -N matlabScript #
matlab -nosplash -nodisplay -r "sbmoptisim_block(150, 3, 0.02, 0.4, 1, 100, 12); quit;"
echo "" 
echo "Done at " `date` 