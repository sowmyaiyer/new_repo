for f in  ../bsubFiles/hisat2_*_S*.bsub
do
	jobname=`basename $f | sed 's/\.bsub//g'`
	cat $HOME/bsub_options_long $f | sed "s/InsertJobName/${jobname}/g" | bsub
done
