#for f in ../bsubFiles/new_patient_samples_*bsub
for f in   ../bsubFiles/chip_retest2_anetta_diagen*bsub
do
	jobname=`basename $f | sed 's/\.bsub//g'`
	cat $HOME/bsub_options_long $f | sed "s/InsertJobName/${jobname}/g" | bsub
done
