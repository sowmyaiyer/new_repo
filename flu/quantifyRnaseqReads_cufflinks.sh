for time in {"6h","12h","24h","48h","mock12h"}
do
        echo """
        module load bowtie2/2-2.1.0
	module load cufflinks/2.2.1
        mkdir  /home/si14w/gnearline/flu/cufflinks_out_${time}
        cufflinks -p 8 -u -g ../txt/hg19_ucsc_genes_and_flu.longer_than_200.gtf --library-type fr-firststrand -o /home/si14w/gnearline/flu/cufflinks_out_${time} /home/si14w/gnearline/flu/tophat_out_${time}/accepted_hits.bam 
        """ > ../bsubFiles/runcufflinks_rnaseq.${time}.bsub
done
