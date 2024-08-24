#!/bin/sh
#PBS -q hotel
#PBS -N run_brian_bk6
#PBS -l nodes=1:ppn=6,walltime=10:00:00
#PBS -o run_brian_bk6.o
#PBS -e run_brian_bk6.e
#PBS -V
#PBS -M shz334@ucsd.edu
#PBS -m abe
dest=$HOME'/Downloads/Polarity-master-2/Zid-Lab/'
input=$dest'Analysis/'
cur=$dest'Analysis/Inputs/'
# Set file names and paths for use in bowtie and data extraction
mapDir='Maps/'
bowDir='bowtie-results/'
# Set index basenames
rdnaDir='Analysis/Indexes/rDNA'
geneDir='Analysis/Indexes/genome'
ecoliDir='Analysis/Indexes/e_coli'
# Set bowtie output filenames
rDNA=$dest$bowDir'rDNA-'
coding=$dest$bowDir'coding-'
map='aligned-'
dchrom=$input'Chrom/'
trimf='Analysis/Inputs/trim/'
sortedGeneDir='SortedFiles/'
trimpy='trimPolyA.py'
fct='FeatureCounts.py'
fctsam='FeatureCounts_sam.py'
# Create bowtie argument filepaths
arg1=$dest$rdnaDir
arg4=$dest$geneDir

# the fastq file contains the sequence for both rDNA and genome.
for f in $input/*.fastq; do
        #To get the name of the file
        file=`basename $f`;
        echo $file;
        #python3 $cur$trimpy:currentdirtrimPolyA.py (to call the function) $dest$trimf$file(param): the path to the fastq file in trim/ $f(param) the file name
        #python3 $cur$trimpy $f $dest$trimf$file ;
        #echo $dest$mapDir$file;
        #all bowtie files in $PATH variable so it can be executed no matter what dir I am in.        
        export PATH=$PATH:$HOME'Downloads/bowtie-master/bowtie' 
        #cutadapt -a CTGTAGGCACCATCAAT -o output.fastq SRR5008134.fastq
        # -a --best make bowtie prints out all the alignments in best-worst-order. --strata makes bowtie only prints out the alignments with the least mismatches.
        # --un <filename> writes all the reads that couldn't be read to the file.
        # $arg1 provides the rDNA indexes. $dest$trimf$file provides the path to the sequence that's trimed by trimPolyA.py
        # $rDNA$file defines the output: output would be in Downloads/.../Zid-Lab/bowtie-results/rDNA-<name of fastq file>
        # --un $coding$file: unmatched reads will be written in Downloads/.../Zid-Lab/bowtie-results/coding-fastqName.
        bowtie -a --best --strata $arg1s $dest'Analysis/Inputs/output.fastq' $rDNA$file --un $coding$file
        # -m <int> suppress all alignments for a particular read if more than <int> alignments exit.
        # $arg4 provides the gene indexes. $coding$file is the file contains all the mismatched reads fro last step.
        # the output map would be in Downloads/.../Zid-Lab/bowtie-results/Maps/<fastq name>
        bowtie -m 1 --best --strata $arg4 $coding$file $dest$mapDir$file
        # call FeatureCounts.py to sort and map all the files under Analysis/chrom.
        # $dest$sortedGeneDir provides .../Zid-Lab/SortedFiles/<fastq name> $dest$mapDir$file provide .../Zid-Lab/Maps/<fastq name>
        # $dest$dchrom provides path to all the files under Analysis/chrom/.
        python3 $cur$fct $dest$sortedGeneDir $dest$mapDir$file $dest$dchrom

done
for f in $input/*.sam; do
	file=`basename $f`;
        echo $file;
	python3 $cur$fctsam $dest$sortedGeneDir $f $dest$dchrom
#done

# do cutadap first to trim off the adapter and then do the two bowtie commands
# may try to delete the .gz files.
# align sort count