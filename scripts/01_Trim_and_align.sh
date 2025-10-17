#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=99G
#SBATCH --export=NONE
#SBATCH --array=4591-6000
#SBATCH --job-name=Trim_and_align_all
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module add Trimmomatic/0.39-Java-11
module add BWA/0.7.17-GCC-10.3.0
module add samblaster/0.1.26-GCC-10.3.0
module add SAMtools/1.13-GCC-10.3.0
module add Sambamba/0.8.2-GCC-10.3.0

echo "start time"
date

#WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_calling

WD=/storage/research/iee_temp_dj20y461/DV_calling

ID=$SLURM_ARRAY_TASK_ID
ITER_FILE=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/sample_metadata/sample_paths_to_realign.txt

SAMPLE_FASTQ_PATH=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITER_FILE)
SAMPLE_NAME=$(echo $SAMPLE_FASTQ_PATH | rev | cut -f1 -d'/' | rev)

######### Trimming ##############################################################################################################################################################################################

ADAPTERS=/software.9/software/Trimmomatic/0.39-Java-11/adapters/NexteraPE-PE.fa

TRIMMED_OUTDIR=$WD/trimmed

if [ ! -d "$TRIMMED_OUTDIR" ]
then
    mkdir $TRIMMED_OUTDIR
fi

## OUT FILES

SAMPLE_R1_trimmed=$TRIMMED_OUTDIR/${SAMPLE_NAME}.R1.trimmed.fastq.gz
SAMPLE_R1_unpaired=$TRIMMED_OUTDIR/${SAMPLE_NAME}.R1.unpaired.fastq.gz
SAMPLE_R2_trimmed=$TRIMMED_OUTDIR/${SAMPLE_NAME}.R2.trimmed.fastq.gz
SAMPLE_R2_unpaired=$TRIMMED_OUTDIR/${SAMPLE_NAME}.R2.unpaired.fastq.gz

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
	-threads 20 \
        $SAMPLE_FASTQ_PATH*_R1_*fastq.gz $SAMPLE_FASTQ_PATH*_R2_*fastq.gz \
        $SAMPLE_R1_trimmed $SAMPLE_R1_unpaired \
        $SAMPLE_R2_trimmed $SAMPLE_R2_unpaired \
        ILLUMINACLIP:${ADAPTERS}:3:30:10 SLIDINGWINDOW:4:20 MINLEN:25 HEADCROP:7


######### Aligning ##############################################################################################################################################################################################

DATA_DIR=$WD/trimmed

GENOME_IDX=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/ref/BWA/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY

## Get the read group info to add to bams from the first read of every file

HEADER=$(zcat $SAMPLE_R1_trimmed | head -n1) ## get first read ID line 

INSTRUMENT=$(echo $HEADER | cut -f1 -d':'| sed 's/@//g')
RUN=$(echo $HEADER | cut -f2 -d':')
FLOWCELL=$(echo $HEADER | cut -f3 -d':')
LANE=$(echo $HEADER | cut -f4 -d':')
BARCODES=$(echo $HEADER | cut -f2 -d' '| cut -f4 -d':' | sed 's/+/-/g')

echo "@RG\tID:${SAMPLE_NAME}\tPL:ILLUMINA\tPU:${FLOWCELL}.${LANE}\tSM:${SAMPLE_NAME}\tCN:MCGILL\tBC:${BARCODES}"

echo $READGROUP_STRING

BAMDIR=$WD/bams/

if [ ! -d "$BAMDIR" ]; then
   mkdir $BAMDIR
fi

BAM_FIX_COORDSORT=${BAMDIR}/${SAMPLE_NAME}.fixmate.coordsorted.bam

bwa mem -t 20 \
        -R $(echo "@RG\tID:${SAMPLE_NAME}\tPL:ILLUMINA\tPU:${FLOWCELL}.${LANE}\tSM:${SAMPLE_NAME}\tCN:MCGILL\tBC:${BARCODES}") \
        $GENOME_IDX \
        $SAMPLE_R1_trimmed \
        $SAMPLE_R2_trimmed | \
samblaster | \
samtools fixmate \
	-m \
        -@ 20 \
        - \
	- | \
samtools sort \
	-@ 20 \
	-O BAM \
	-o $BAM_FIX_COORDSORT \
	-

samtools index -c -@ 20 $BAM_FIX_COORDSORT

## End time
echo "end time"
date

##################################################################
####### calculate some stats ####################################
#################################################################

samtools stats $BAM_FIX_COORDSORT > ${BAMDIR}/${SAMPLE_NAME}.stats
samtools coverage $BAM_FIX_COORDSORT  > ${BAMDIR}/${SAMPLE_NAME}.depths

## Calculate some averages and add to the end of each file
grep 'chromosome' ${BAMDIR}/${SAMPLE_NAME}.depths | awk '{ sum += $6; n++ } END { if (n > 0) print "\n## Avg. perc. bases. covered (chroms) = " sum / n; }' >> ${BAMDIR}/${SAMPLE_NAME}.depths
grep 'chromosome' ${BAMDIR}/${SAMPLE_NAME}.depths | awk '{ sum += $7; n++ } END { if (n > 0) print "## Avg. coverage (chroms) = " sum / n; }' >> ${BAMDIR}/${SAMPLE_NAME}.depths


######### Cleanup ##############################################################################################################################################################################################

## the raw bam files probably don't need to be kept. . . . lets see how it goes for now. 

## Can I store bams as CRAM? 



