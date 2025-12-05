#!/bin/bash
#SBATCH --job-name=jasmine_filter
#SBATCH --output=/storage/homefs/kw23y068/logfiles/jasmine_filter_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/jasmine_filter_%j.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --mem=64G
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

PROJECT_DIR=/storage/scratch/iee_evol/kw23y068
OUTPUT_DIR=${PROJECT_DIR}/jasmine_${START_PAD}_${END_PAD}

module -q reset
module load HTSlib/1.12-GCC-10.3.0
module load BCFtools/1.12-GCC-10.3.0
module load VCFtools/0.1.16-GCC-10.3.0

BASE="Manta_Smoove_Delly"
VCFlist=${OUTPUT_DIR}/${BASE}.list
OUTfile=${OUTPUT_DIR}/${BASE}.vcf
STORAGE_DIR="/storage/research/iee_evol/Korbi"
REFERENCE="${STORAGE_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna"

cd ${OUTPUT_DIR}

bcftools query -f '%CHROM %POS %INFO/END %INFO/SVTYPE %INFO/SVLEN %INFO/SUPP_VEC %INFO/SUPP\n' ${OUTfile} > ${BASE}.bed

echo "nb of SV"
grep -v "#" ${BASE}.vcf | wc -l 

perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/' ${BASE}.vcf | sed -e 's/\(.\)/\1 /g' > ${BASE}.overlap.txt

#remove the ones on Unplaced (no validity) and W (only covered on half) and MT and Z
# vcftools --vcf ${BASE}.vcf \
# --not-chr Chr_W \
# --not-chr Chr_Unplaced  \
# --not-chr Chr_MT \
# --not-chr Chr_Z \
# --remove-filtered-all --recode --recode-INFO-all --stdout | gzip > ${BASE}.Autosomes.vcf.gz