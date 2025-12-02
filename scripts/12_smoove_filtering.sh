#!/bin/bash
#SBATCH --job-name=smoove_filter
#SBATCH --output=/storage/homefs/kw23y068/logfiles/smoove_filter_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/smoove_filter_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=epyc2
#SBATCH --account=gratis
#SBATCH --wckey=noop

####################
# PARAMETERS TO SET
####################
# May be altered for feeding in arguments to script
# Define which batch of samples to process
BATCH_START=1    # Change to 81, 161, 241, 321, 401 for other batches
BATCH_END=80     # Change to 160, 240, 320, 400, 480 for other batches

START_PAD=$(printf "%03d" "${BATCH_START}")
END_PAD=$(printf "%03d" "${BATCH_END}")

echo "Start smoove filtering"
echo "Date: $(date)"

# Define paths
PROJECT_DIR="/storage/research/iee_evol/Korbi"
REFERENCE="${PROJECT_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna"
OUTPUT_DIR="/storage/scratch/iee_evol/kw23y068/smoove_${START_PAD}_${END_PAD}"
#/storage/scratch/iee_evol/kw23y068/smoove_001_080/smoove_merged_genotyped_001_080-smoove.genotyped.vcf.gz
#/storage/scratch/iee_evol/kw23y068/smoove_001_080/smoove_merged_genotyped_001_080-smoove.genotyped.vcf
MERGED_VCF="${OUTPUT_DIR}/smoove_merged_genotyped_${START_PAD}_${END_PAD}-smoove.genotyped.vcf"
FILTERED_VCF="${OUTPUT_DIR}/smoove_${START_PAD}_${END_PAD}_merged.filtered.vcf.gz"
FILTERED_1_VCF="${OUTPUT_DIR}/smoove_${START_PAD}_${END_PAD}_merged.filtered1.vcf.gz"
FILTERED_1_SORTED_VCF="${OUTPUT_DIR}/smoove_${START_PAD}_${END_PAD}_merged.filter1.sorted.vcf.gz"

# Check if input file exists
if [ ! -f "${MERGED_VCF}" ]; then
    echo "ERROR: Input file does not exist: ${MERGED_VCF}"
    echo "Check if previous scripts ran correctly."
    exit 1
fi

cp ${MERGED_VCF} ${FILTERED_VCF}

echo "nb of SV detected by Smoove"
grep -v "#" ${FILTERED_VCF} | wc -l 

# Load required modules (we use slightly newer)
module -q reset
# module load HTSlib/1.9-intel-2018b
# module load BCFtools/1.9-intel-2018b
# module load VCFtools/0.1.16-intel-2018b-Perl-5.28.0
module load HTSlib/1.12-GCC-10.3.0
module load BCFtools/1.12-GCC-10.3.0
module load VCFtools/0.1.16-GCC-10.3.0

bgzip -c ${FILTERED_VCF}

#Filter by developers
bcftools view -i '(SVTYPE = "DEL" & FMT/DHFFC[0] < 0.7) | (SVTYPE = "DUP" & FMT/DHBFC[0] > 1.3) | (SVTYPE = "INV") | (SVTYPE = "INS")' \
    -O z \
    -o ${FILTERED_1_VCF} \
    ${FILTERED_VCF}
bcftools sort \
    -O z \
    -o ${FILTERED_1_SORTED_VCF} \
    ${FILTERED_1_VCF} ;

zcat ${FILTERED_1_SORTED_VCF} | grep -v "#" | wc -l 

# We overwrite the Variables of the authors with ours
INPUT_VCF=${FILTERED_1_SORTED_VCF} # we will try not to modify this one
VCF_FOLDER=${OUTPUT_DIR}/TMP_VCF
mkdir -p ${VCF_FOLDER}
OUTPUT_VCF="${OUTPUT_DIR}/smoove_${START_PAD}_${END_PAD}_filtered.vcf"
RScripts="/storage/homefs/kw23y068/software/scripts/R_smoove"
RefFasta=${REFERENCE}

cp $INPUT_VCF ${VCF_FOLDER}/raw.vcf.gz
gunzip ${VCF_FOLDER}/raw.vcf.gz

#filter vcf -i (include, -O vcf format -o
bcftools filter -i'INFO/SVLEN>=50 | INFO/SVLEN<=-50' \
    -o ${VCF_FOLDER}/raw_sorted.noTRA_50bp.vcf \
    -Ov ${VCF_FOLDER}/raw.vcf
echo "total number of SVs > 50b"
grep -v "#" ${VCF_FOLDER}/raw_sorted.noTRA_50bp.vcf \
    | wc -l

# All thanks to Marc-André Lemay and the people I copy from second hand
#then we use the reference to get the sequence with a R scripts graciously provided by Marc-André Lemay
module -q reset
module load Anaconda3/2024.02-1
source $(conda info --base)/etc/profile.d/conda.sh
conda activate R_all
Rscript ${RScripts}/add_explicit_seq.r ${VCF_FOLDER}/raw_sorted.noTRA_50bp.vcf ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ.vcf ${RefFasta}

conda deactivate

module -q reset
module load HTSlib/1.12-GCC-10.3.0
module load BCFtools/1.12-GCC-10.3.0
module load VCFtools/0.1.16-GCC-10.3.0
#Export sequences for advanced filtering
bcftools query \
    -f '%CHROM %POS %INFO/END %INFO/SVTYPE %INFO/SVLEN %REF %ALT\n' \
    ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ.vcf \
    > ${VCF_FOLDER}/SV_data_with_seq.txt

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
bgzip -c ${VCF_FOLDER}/blacklist.bed \
    > ${VCF_FOLDER}/blacklist.bed.gz
tabix -s1 -b2 -e2 ${VCF_FOLDER}/blacklist.bed.gz

#remove blacklist of variants
bcftools view -T \
    ^${VCF_FOLDER}/blacklist.bed.gz \
    ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ.vcf \
    > ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ_Nfiltered.vcf
echo "SVs after filtration for N seq" 
grep -v "#" \
    ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ_Nfiltered.vcf \
    | wc -l 

#keep the filtered vcf
cp ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ_Nfiltered.vcf ${OUTPUT_VCF}
bgzip -c ${OUTPUT_VCF} > ${OUTPUT_VCF}.gz
tabix ${OUTPUT_VCF}.gz

echo "End smoove filtering"
echo "Date: $(date)"
