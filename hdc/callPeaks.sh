bamFile_treat=$1
bamFile_control=$2
output_name=$3
macs2 callpeak -t ${bamFile_treat}  -c ${bamFile_control} -B --SPMR -g hs --nomodel --extsize 200 -n $HOME/gnearline/hdc/macs2_out/${output_name}
#macs2 bdgcmp -t  $HOME/gnearline/hdc/macs2_out/${output_name}_treat_pileup.bdg -c  $HOME/gnearline/hdc/macs2_out/${output_name}_control_lambda.bdg -o  $HOME/gnearline/hdc/macs2_out/${output_name}_single.bdg -m FE
#$HOME/gnearline/common_scripts/bdg2bw  $HOME/gnearline/hdc/macs2_out/${output_name}_single.bdg $HOME/gnearline/common_scripts/hg19.genome
