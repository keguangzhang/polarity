export PATH=$PATH:$HOME'path/to/cutadapt'
cutadapt -a CTGTAGGCACCATCAAT -m 18 -M 32 -o Inputs/aftercutadapt/outputXXX.fastq Inputs/inputfastq/SRRXXXXXX.fastq.gz
# Cutadapt cut adaptor sequence. See cutadapt manual for command line and option description.

export PATH=$PATH:$HOME'path/to/bowtie' 
bowtie -y -a --best --strata Indexes/rDNA Inputs/aftercutadapt/outputXXX.fastq Inputs/uselessXXX.fastq --un Inputs/afterbowtie1/XXX.fastq 
bowtie --chunkmbs 512 -y -a -m 1 -p 8 --best --strata Indexes/genome Inputs/afterbowtie1/XXX.fastq Inputs/afterbowtie2/XXX.fastq
# Align sequence to rDNA and genome respectively not reporting reads that have more then 1 alignments. 
# See cutadapt manual for command line and option description

python3 Inputs/FeatureCounts.py Inputs/afterfeaturecounts/XXX Inputs/afterbowtie2/XXX.fastq Inputs/Chrom/
#align sequence to genes and count number

python3 Inputs/polarity.py Inputs/afterfeaturecounts/ Inputs/afterpolarity/
#calculate polarity