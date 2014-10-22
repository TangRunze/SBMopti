# 
#$ -cwd 
#$ -j y 
#$ -pe orte 12
#$ -S /bin/bash 
#
# Name the job #$ -N matlabScript #
matlab -nosplash -nodisplay -r "sbmoptisim(300, 3, 0.1, 1, 1000, 12); quit;"
echo "" 
echo "Done at " `date` 