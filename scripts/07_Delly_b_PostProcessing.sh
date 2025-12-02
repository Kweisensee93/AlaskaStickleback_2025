#!/bin/bash
#SBATCH --job-name=Delly_PostProcessing
#SBATCH --output=/storage/homefs/kw23y068/logfiles/Delly_PostProcessing_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/Delly_PostProcessing_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=epyc2

####################
# PARAMETERS TO SET
####################
# May be altered for feeding in arguments to script
# Define which batch of samples to process
BATCH_START=1    # Change to 81, 161, 241, 321, 401 for other batches
BATCH_END=80     # Change to 160, 240, 320, 400, 480 for other batches

# Handle with care: bgzip and tabix are used with -f to overwrite existing files
# We use module reset instead of module purge on our cluster.
# If you want to supress the warnings that default modules are empty, update script
# module reset --> module -q reset

## We use newer versions available on cluster
#module load HTSlib/1.9-intel-2018b
#module load BCFtools/1.9-intel-2018b
#module load VCFtools/0.1.16-intel-2018b-Perl-5.28.0
module load HTSlib/1.12-GCC-10.3.0
module load BCFtools/1.12-GCC-10.3.0
module load VCFtools/0.1.16-GCC-10.3.0

# we already have bgzip from 07_Delly_a script
#bgzip SPI_BJO_ROS_sites.vcf
# bgzip ${OUTPUT_VCF}
# tabix -p vcf ${OUTPUT_VCF}.gz

PROJECT_DIR=/storage/research/iee_evol/Korbi
OUTPUT_DIR=/storage/scratch/iee_evol/kw23y068/

START_PAD=$(printf "%03d" "${BATCH_START}")
END_PAD=$(printf "%03d" "${BATCH_END}")
MERGE_DIR="${OUTPUT_DIR}/Delly_Merge_${START_PAD}_${END_PAD}"

#INPUT_VCF=SPI_BJO_ROS_sites.vcf.gz # we will try not to modify this one
INPUT_VCF=${MERGE_DIR}/Delly_merged_${START_PAD}_${END_PAD}.vcf.gz
VCF_FOLDER=${MERGE_DIR}/TMP_VCF
OUTPUT_VCF=${MERGE_DIR}/Delly_merged_filtered_${START_PAD}_${END_PAD}.vcf

# we need to adapt the R scripts
#RScripts=SV_Merot/Rscripts
# I adapted the following 3 R scripts from
## https://github.com/clairemerot/SR_SV/blob/main/01_scripts/Rscripts/
# add_explicit_seq_delly.r
# add_info_bcf.r
# fix_sniffles_delly.r
RScripts=/storage/homefs/kw23y068/software/scripts/R_Delly
#RefFasta=SV_Merot/Ref/Puffin_ReferenceGenome.NU.MT.fasta
# I only have .fna available --> The script ran fine with .fna --> no need to convert
REFERENCE=/storage/research/iee_evol/Korbi/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna

# Check if R scripts and reference exist
if [ ! -d "${RScripts}" ]; then
    echo "ERROR: R scripts directory not found: ${RScripts}"
    exit 1
fi

if [ ! -f "${REFERENCE}" ]; then
    echo "ERROR: Reference fasta not found: ${REFERENCE}"
    exit 1
fi

if [ ! -f "${INPUT_VCF}" ]; then
    echo "ERROR: Input VCF not found: ${INPUT_VCF}"
    exit 1
fi

mkdir -p ${VCF_FOLDER}
cp $INPUT_VCF ${VCF_FOLDER}/raw.vcf.gz
gunzip ${VCF_FOLDER}/raw.vcf.gz

#check what it look like
#grep -v ^\#\# $VCF_FOLDER/raw.vcf | head 
#tail $VCF_FOLDER/raw.vcf 
echo "total number of SVs"
grep -v "#" ${VCF_FOLDER}/raw.vcf | wc -l  

#filter out TRA & BND & INS
bcftools filter -i'INFO/SVTYPE!="TRA" & INFO/SVTYPE!="BND" & INFO/SVTYPE!="INS"' \
    -o ${VCF_FOLDER}/raw.noTRA.vcf -Ov ${VCF_FOLDER}/raw.vcf 
#grep -v ^\#\# $VCF_FOLDER/raw_sorted.noTRA.vcf | head
echo "total number of SVs restricted to  DEL, INV, DUP"
grep -v "#" ${VCF_FOLDER}/raw.noTRA.vcf | wc -l 

#Keep INS and add the sequence of INS
bcftools filter -i 'INFO/SVTYPE=="INS"' \
    -o ${VCF_FOLDER}/raw.INS.vcf -Ov ${VCF_FOLDER}/raw.vcf 
echo "total number of INS SVs"
grep -v "#" ${VCF_FOLDER}/raw.INS.vcf | wc -l 

#put the field consensus in ALT
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%INFO/CONSENSUS\t%QUAL\t%FILTER\n' \
    ${VCF_FOLDER}/raw.INS.vcf > ${VCF_FOLDER}/raw.INS.info1
cat ${VCF_FOLDER}/raw.INS.vcf | grep -v "#" | cut -f8 > ${VCF_FOLDER}/raw.INS.info2
paste ${VCF_FOLDER}/raw.INS.info1 ${VCF_FOLDER}/raw.INS.info2 > ${VCF_FOLDER}/raw.INS.info

#re-merge everything
echo "Merging INS back with other SV types..."
grep ^"#" ${VCF_FOLDER}/raw.vcf > ${VCF_FOLDER}/raw_sorted.noTRA.vcf
(grep -v ^"#" ${VCF_FOLDER}/raw.noTRA.vcf; grep -v ^"#" ${VCF_FOLDER}/raw.INS.info) | \
    sort -k1,1 -k2,2n >> ${VCF_FOLDER}/raw_sorted.noTRA.vcf

#grep -v ^\#\# $VCF_FOLDER/raw_sorted.vcf | head
echo "total number of SVs"
grep -v "#" ${VCF_FOLDER}/raw_sorted.noTRA.vcf | wc -l 

#we need to add a field for SVLEN
##step 1 export position 
bcftools query -f '%CHROM\t%POS\t%INFO/END\n' \
    ${VCF_FOLDER}/raw_sorted.noTRA.vcf > ${VCF_FOLDER}/raw_sorted.noTRA.info

#step 2 calculate length
# On cluster we use reset instead of purge
# We set up R in Conda environment instead of module (works better - TRUST!)
#module purge
#module load R/4.1.0-foss-2021a

module -q reset
module load Anaconda3/2024.02-1
source $(conda info --base)/etc/profile.d/conda.sh
conda activate R_all

Rscript ${RScripts}/add_info_bcf.r ${VCF_FOLDER}/raw_sorted.noTRA.info 

conda deactivate
module -q reset

module load HTSlib/1.12-GCC-10.3.0
module load BCFtools/1.12-GCC-10.3.0
module load VCFtools/0.1.16-GCC-10.3.0

bgzip -f ${VCF_FOLDER}/raw_sorted.noTRA.info.annot
tabix -f -s1 -b2 -e2 ${VCF_FOLDER}/raw_sorted.noTRA.info.annot.gz

##step3 prepare the header
echo -e '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">' \
    > ${VCF_FOLDER}/raw_sorted.noTRA.info.annot.hdr

##step4 run bcftools annotate
#-a is the annotation file (tabix and bgzip, it needs at least CHROM and POS, -h are the header lines to add, -c are the meaning of the column in the annotation file
bcftools annotate -a ${VCF_FOLDER}/raw_sorted.noTRA.info.annot.gz \
    -h ${VCF_FOLDER}/raw_sorted.noTRA.info.annot.hdr \
    -c CHROM,POS,INFO/SVLEN \
    ${VCF_FOLDER}/raw_sorted.noTRA.vcf \
    > ${VCF_FOLDER}/raw_sorted.noTRA_bis.vcf

#filter vcf -i (include, -O vcf format -o
echo "Filtering SVs by size (≥50bp)..."
bcftools filter -i'INFO/SVLEN>=50 | INFO/SVLEN<=-50' \
    -o ${VCF_FOLDER}/raw_sorted.noTRA_50bp.vcf \
    -Ov ${VCF_FOLDER}/raw_sorted.noTRA_bis.vcf
echo "total number of SVs > 50b"
grep -v ^\#\# ${VCF_FOLDER}/raw_sorted.noTRA_50bp.vcf | wc -l 

#then we use the reference to get the sequence with a R scripts graciously provided by Marc-André Lemay
module reset
module load Anaconda3/2024.02-1
source $(conda info --base)/etc/profile.d/conda.sh
conda activate R_all

Rscript ${RScripts}/add_explicit_seq_delly.r \
    ${VCF_FOLDER}/raw_sorted.noTRA_50bp.vcf \
    ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ.vcf \
    ${REFERENCE}

conda deactivate
module reset
#Export sequences for advanced filtering

module load HTSlib/1.12-GCC-10.3.0
module load BCFtools/1.12-GCC-10.3.0
module load VCFtools/0.1.16-GCC-10.3.0

bcftools query -f '%CHROM %POS %INFO/END %INFO/SVTYPE %INFO/SVLEN %REF %ALT\n' \
    ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ.vcf \
    > ${VCF_FOLDER}/SV_data_with_seq.txt

# Create blacklists
echo "Creating blacklists..."

#blacklist because of N string > 10 (possible junction of contigs 
grep -P "N{10,}" ${VCF_FOLDER}/SV_data_with_seq.txt \
    | awk '{print $1 "\t" $2 "\t" $6 "\t" $7}' \
    > ${VCF_FOLDER}/N10_blacklist.bed
echo "SVs excluded because of >10N" 
wc -l ${VCF_FOLDER}/N10_blacklist.bed 

#blacklist because missing seq
cat ${VCF_FOLDER}/SV_data_with_seq.txt \
    | awk '{if ($6 == "N") print $1 "\t" $2 "\t" $6 "\t" $7;}' \
    > ${VCF_FOLDER}/N_blacklist.bed
echo "SVs excluded because absence of sequence ref" 
wc -l ${VCF_FOLDER}/N_blacklist.bed 

#blacklist because missing seq
cat  ${VCF_FOLDER}/SV_data_with_seq.txt \
    | awk '{print $1 "\t" $2 "\t" $6 "\t" $7}' \
    | grep -P "<" > ${VCF_FOLDER}/N_blacklist_bis.bed
echo "SVs excluded because absence of sequence alt" 
wc -l ${VCF_FOLDER}/N_blacklist_bis.bed 

#full blacklist
cat ${VCF_FOLDER}/N_blacklist.bed \
    ${VCF_FOLDER}/N_blacklist_bis.bed \
    ${VCF_FOLDER}/N10_blacklist.bed \
    | sort -k1,1 -k2,2n \
    > ${VCF_FOLDER}/blacklist.bed
#head $VCF_FOLDER/blacklist.bed
bgzip -f -c ${VCF_FOLDER}/blacklist.bed > ${VCF_FOLDER}/blacklist.bed.gz
tabix -f -s1 -b2 -e2 ${VCF_FOLDER}/blacklist.bed.gz

#remove blacklist of variants
bcftools view -T ^${VCF_FOLDER}/blacklist.bed.gz \
    ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ.vcf \
    > ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ_Nfiltered.vcf
echo "SVs after filtration for N seq" 
grep -v "#" ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ_Nfiltered.vcf | wc -l 

#keep the filtered vcf
cp ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ_Nfiltered.vcf ${OUTPUT_VCF}
bgzip -f -c ${OUTPUT_VCF} > ${OUTPUT_VCF}.gz
tabix -f ${OUTPUT_VCF}.gz
