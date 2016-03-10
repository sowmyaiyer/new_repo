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

echo """
module load python/2.7.9
module load blas/08_2013
module load atlas/3.10.2
module load mpich/3.0.4
module load gcc/5.1.0
module load lapack/3.5.0_gcc_5.1.0
module load intel/mkl_libraries_from_composer_xe_2013_sp1
cd ~/NucleoATAC
source $HOME/venv/bin/activate
nucleoatac run --bed /home/si14w/gnearline/hdc/HumanDC_ATAC_broad/ATAC_broad_merged_E70_E72_4h_LPS_peaks.broadPeak --bam /project/umw_garberlab/human/DC/ATAC-Seq/ATAC_ChIP_12_22_15/pipe_E70_E72_merged/alignments/merged_E70_E72_4h_LPS.cutadapt.sorted.no_dups.filt.bam --fasta /share/data/umw_biocore/Genomes/human/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa --out /home/si14w/gnearline/hdc/nucleoatac_out/nucleoatac_out_E70_E72_4h_LPS
""" > ../bsubFiles/nucleoatac_merged_E70_E72_4h_LPS.bsub

echo """
module load python/2.7.9
module load blas/08_2013
module load atlas/3.10.2
module load mpich/3.0.4
module load gcc/5.1.0
module load lapack/3.5.0_gcc_5.1.0
module load intel/mkl_libraries_from_composer_xe_2013_sp1
cd ~/NucleoATAC
source $HOME/venv/bin/activate
nucleoatac run --bed /home/si14w/gnearline/hdc/HumanDC_ATAC_broad/ATAC_broad_merged_E70_E72_2h_LPS_peaks.broadPeak --bam /project/umw_garberlab/human/DC/ATAC-Seq/ATAC_ChIP_12_22_15/pipe_E70_E72_merged/alignments/merged_E70_E72_2h_LPS.cutadapt.sorted.no_dups.filt.bam --fasta /share/data/umw_biocore/Genomes/human/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa --out /home/si14w/gnearline/hdc/nucleoatac_out/nucleoatac_out_E70_E72_2h_LPS
""" > ../bsubFiles/nucleoatac_merged_E70_E72_2h_LPS.bsub


echo """
module load python/2.7.9
module load blas/08_2013
module load atlas/3.10.2
module load mpich/3.0.4
module load gcc/5.1.0
module load lapack/3.5.0_gcc_5.1.0
module load intel/mkl_libraries_from_composer_xe_2013_sp1
cd ~/NucleoATAC
source $HOME/venv/bin/activate
nucleoatac run --bed /home/si14w/gnearline/hdc/HumanDC_ATAC_broad/ATAC_broad_merged_E70_E72_30min_LPS_peaks.broadPeak --bam /project/umw_garberlab/human/DC/ATAC-Seq/ATAC_ChIP_12_22_15/pipe_E70_E72_merged/alignments/merged_E70_E72_30min_LPS.cutadapt.sorted.no_dups.filt.bam --fasta /share/data/umw_biocore/Genomes/human/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa --out /home/si14w/gnearline/hdc/nucleoatac_out/nucleoatac_out_E70_E72_30min_LPS
""" > ../bsubFiles/nucleoatac_merged_E70_E72_30min_LPS.bsub

echo """
module load python/2.7.9
module load blas/08_2013
module load atlas/3.10.2
module load mpich/3.0.4
module load gcc/5.1.0
module load lapack/3.5.0_gcc_5.1.0
module load intel/mkl_libraries_from_composer_xe_2013_sp1
cd ~/NucleoATAC
source $HOME/venv/bin/activate
nucleoatac run --bed /home/si14w/gnearline/hdc/HumanDC_ATAC_broad/ATAC_broad_merged_E70_E72_0h_LPS_peaks.broadPeak --bam /project/umw_garberlab/human/DC/ATAC-Seq/ATAC_ChIP_12_22_15/pipe_E70_E72_merged/alignments/merged_E70_E72_0h_LPS.cutadapt.sorted.no_dups.filt.bam --fasta /share/data/umw_biocore/Genomes/human/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa --out /home/si14w/gnearline/hdc/nucleoatac_out/nucleoatac_out_E70_E72_0h_LPS
""" > ../bsubFiles/nucleoatac_merged_E70_E72_0h_LPS.bsub


echo """
module load python/2.7.9
module load blas/08_2013
module load atlas/3.10.2
module load mpich/3.0.4
module load gcc/5.1.0
module load lapack/3.5.0_gcc_5.1.0
module load intel/mkl_libraries_from_composer_xe_2013_sp1
cd ~/NucleoATAC
source $HOME/venv/bin/activate
nucleoatac run --bed /home/si14w/gnearline/hdc/HumanDC_ATAC_broad/ATAC_broad_merged_E70_E72_24h_LPS_peaks.broadPeak --bam /project/umw_garberlab/human/DC/ATAC-Seq/ATAC_ChIP_12_22_15/pipe_E70_E72_merged/alignments/merged_E70_E72_24h_LPS.cutadapt.sorted.no_dups.filt.bam --fasta /share/data/umw_biocore/Genomes/human/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa --out /home/si14w/gnearline/hdc/nucleoatac_out/nucleoatac_out_E70_E72_24h_LPS
""" > ../bsubFiles/nucleoatac_merged_E70_E72_24h_LPS.bsub
