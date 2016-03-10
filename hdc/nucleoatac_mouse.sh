# one time
virtualenv /home/si14w/venv/
module load python/2.7.9
cd /
source $HOME/venv/bin/activate
# Now you are in virtual environment

module load blas/08_2013
module load atlas/3.10.2
module load mpich/3.0.4
module load gcc/5.1.0
module load lapack/3.5.0_gcc_5.1.0
module load intel/mkl_libraries_from_composer_xe_2013_sp1
pip install numpy --upgrade --force -vv
git clone https://github.com/GreenleafLab/NucleoATAC.git
cd NucleoATAC
pip install . --upgrade

# At one time, even after loading atlas/3.10.2, numpy installation was pointing to a .so file in atlas/3.1.0. Error was "ImportError: /share/pkg/atlas/3.10.1/lib/libtatlas.so: undefined symbol: clapack_ilaenv". This error would show up at runtime
# Upon Arjan's suggestion, I did rm -Rf ~/.cache, loaded atlas/3.10.2 again and did pip install numpy --upgrade --force -vv. Took care of it.
for time in {"0h","30m_LPS","1h_LPS","2h_LPS","4h_LPS","6h_LPS"}
do
echo """
#module load python/2.7.9_packages/macs2/2.1.0
#macs2 callpeak --broad -t /project/umw_garberlab/mouse/DC/ATAC-Seq/merge.ATAC_12_15_15_m2/alignments/MM1_MM2_${time}.cutadapt.sorted.no_dups.filt.bam -f BAM -g mm -n mouse_ATAC_broad_M1_M2_${time} --outdir ../MouseDC_ATAC_broad/


module load python/2.7.9
module load blas/08_2013
module load atlas/3.10.2
module load mpich/3.0.4
module load gcc/5.1.0
module load lapack/3.5.0_gcc_5.1.0
module load intel/mkl_libraries_from_composer_xe_2013_sp1

cd ~/NucleoATAC
source $HOME/venv/bin/activate
nucleoatac run --bed /home/si14w/gnearline/hdc/MouseDC_ATAC_broad/mouse_ATAC_broad_M1_M2_${time}_peaks.broadPeak --bam /project/umw_garberlab/mouse/DC/ATAC-Seq/merge.ATAC_12_15_15_m2/alignments/MM1_MM2_${time}.cutadapt.sorted.no_dups.filt.bam --fasta /share/data/umw_biocore/genome_data/mouse/mm10/mm10.fa --out /home/si14w/gnearline/hdc/mouse_nucleoatac_out/nucleoatac_out_M1_M2_${time}
""" > ../bsubFiles/mouse_nucleoatac_merged_MM1_MM2_${time}.bsub
done
