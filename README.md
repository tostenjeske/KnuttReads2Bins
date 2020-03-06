# KnuttReads2Bins

This is the version of the pipeline as it was used for the paper. Some data outputs weren't used and have been disabled. 

## Setup

The memory requirement depends on the selected Kaiju database in the `config.yml` file. If you don't have at least 170 GB available, replace "nr_euk" with "refseq" (44GB).

1. Install the latest version of a conda distribution, sra-tools (SRA Toolkit) and Snakemake (>=5.10)
2. Clone this repository (paper branch)
3. Download the read data into the `input` folder
4. Run the `paper` rule with conda enabled

``` sh
conda create -n snake -c bioconda -c conda-forge snakemake>=5.10 sra-tools
conda activate snake
git clone -b paper https://github.com/KnuttPipeline/KnuttReads2Bins.git
cd KnuttReads2Bins
download_sra(){
  fastq-dump --split-files -O input $1 --gzip
  mv input/$1_1.fastq.gz input/$2_R1.fastq.gz
  mv input/$1_2.fastq.gz input/$2_R2.fastq.gz
}
download_sra SRR11128853 ENR-Ac244days_meta
download_sra SRR11128854 ENR-Ac211days_meta
download_sra SRR11128855 PFL9_meta
snakemake -prj 16 --use-conda paper | tee run.log
# Run with 16 cores available
```
