#!/bin/bash
#SBATCH --job-name=jasmine
#SBATCH --output=/storage/homefs/kw23y068/logfiles/jasmine_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/jasmine_%j.err
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
mkdir -p ${OUTPUT_DIR}
STORAGE_DIR="/storage/research/iee_evol/Korbi"
REFERENCE="${STORAGE_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna"

# Get sample: Go to ouptut directory and copy from SV callers
cd ${OUTPUT_DIR}
# Define expected files
DELLY="${OUTPUT_DIR}/Delly_Merge_${START_PAD}_${END_PAD}/Delly_merged_filtered_${START_PAD}_${END_PAD}.vcf.gz"
MANTA="${OUTPUT_DIR}/manta_${START_PAD}_${END_PAD}/manta_${START_PAD}_${END_PAD}_merged.vcf.gz"
SMOOVE="${OUTPUT_DIR}/smoove_${START_PAD}_${END_PAD}/smoove_${START_PAD}_${END_PAD}_filtered.vcf.gz"

# Safety check - we all know typos and stuff
for f in "$DELLY" "$MANTA" "$SMOOVE" "$REFERENCE"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Required file not found: $f" >&2
        echo "Check pathway and previous callers."
        exit 1
    fi
done

# Copy files - Note: existing files will be overwritten
cp -f "$DELLY" "$MANTA" "$SMOOVE" .

gunzip -f *.vcf.gz

ls *.vcf > Manta_Smoove_Delly.list

# Start with Jasmine

module load Anaconda3/2024.02-1
source $(conda info --base)/etc/profile.d/conda.sh

# I leave out the module reset to avoid issues with conda
#module -q reset
# I leave the module load for BCFtools and SAMtools - may be needed by jasmine ?
#module load BCFtools/1.9-intel-2018b
#module load SAMtools/1.9-GCC-7.3.0-2.30
module load BCFtools/1.12-GCC-10.3.0
module load SAMtools/1.13-GCC-10.3.0

#input a list of vcf files (should not be zip)
BASE="Manta_Smoove_Delly"

VCFlist=${BASE}.list
OUTfile=${BASE}.vcf

# max_dist_linear may be set to 0.1 ?
# Why use min_dist=50 ? to avoid merging very close calls ? --> for smaller SVs, we may need to adapt
# original creators (the ones who got copied by the ones we copy) didnt have -- allow_intrasample
# maybe add --output_genotypes
# reference used --default_zero_genotype \
# reference used --mutual_distance \
jasmine file_list=${VCFlist} \
    out_file=${OUTfile} \
    genome_file=${REFERENCE} \
    --ignore_strand \
    --max_dist_linear=0.5 \
    --max_dist=1000 \
    --min_dist=50 \
    --use_end \
    --normalize_type \
    --allow_intrasample \
    --output_genotypes \
    threads=${SLURM_CPUS_PER_TASK}

conda deactivate
