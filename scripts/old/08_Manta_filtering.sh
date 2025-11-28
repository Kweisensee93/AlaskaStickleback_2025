#!/bin/bash
#SBATCH --job-name=manta_filter
#SBATCH --array=1-80
#SBATCH --output=/storage/homefs/kw23y068/logfiles/manta_filter_%A_%a.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/manta_filter_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=epyc2

# Define paths
SAMPLE_LIST=/storage/research/iee_evol/Korbi/output/bams_full.csv
SAMPLE_PATH=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${SAMPLE_LIST}")
SAMPLE_FILE=$(basename "${SAMPLE_PATH}")
SAMPLE_NAME="${SAMPLE_FILE%.fixmate.coordsorted.bam}"

PROJECT_DIR=/storage/scratch/iee_evol/kw23y068
MANTA_RUN_DIR=${PROJECT_DIR}/Manta/${SAMPLE_NAME}
STORAGE_DIR=/storage/research/iee_evol/Korbi
REFERENCE=${STORAGE_DIR}/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna

# Input VCF from Manta
INPUT_VCF=${MANTA_RUN_DIR}/results/variants/diploidSV.vcf.gz

# Output directory
OUTPUT_DIR=${PROJECT_DIR}/Manta_filtered/${SAMPLE_NAME}
mkdir -p ${OUTPUT_DIR}

VCF_FOLDER=${OUTPUT_DIR}/TMP_VCF
OUTPUT_VCF=${OUTPUT_DIR}/${SAMPLE_NAME}_manta_filtered.vcf

mkdir -p ${VCF_FOLDER}

# Decompress initial VCF
gunzip -c ${INPUT_VCF} > ${VCF_FOLDER}/raw.vcf

echo "Sample: ${SAMPLE_NAME}"
echo "Total number of SVs detected by Manta:"
grep -v "#" ${VCF_FOLDER}/raw.vcf | wc -l

# Load required modules (we use slightly newer)
module -q reset
# module load HTSlib/1.9-intel-2018b
# module load BCFtools/1.9-intel-2018b
# module load VCFtools/0.1.16-intel-2018b-Perl-5.28.0
module load HTSlib/1.12-GCC-10.3.0
module load BCFtools/1.12-GCC-10.3.0
module load VCFtools/0.1.16-GCC-10.3.0

# Convert inversions using Manta's convertInversion.py
cp ${VCF_FOLDER}/raw.vcf ${VCF_FOLDER}/raw.initial.vcf

bgzip ${VCF_FOLDER}/raw.initial.vcf

#module purge
#module load Miniconda3/4.9.2
#conda activate miniconda3_envs/manta

# Path to convertInversion.py - adjust if needed
MANTA_SHARE=/storage/homefs/kw23y068/software/manta-1.6.0.centos6_x86_64/share/manta-1.6.0-1
if [ -f "${MANTA_SHARE}/libexec/convertInversion.py" ]; then
    python ${MANTA_SHARE}/libexec/convertInversion.py \
        ${MANTA_SHARE}/libexec/samtools \
        ${REFERENCE} \
        ${VCF_FOLDER}/raw.initial.vcf > ${VCF_FOLDER}/raw.initial2.vcf
else
    echo "Warning: convertInversion.py not found, skipping inversion conversion"
    cp ${VCF_FOLDER}/raw.initial.vcf ${VCF_FOLDER}/raw.initial2.vcf
fi

conda deactivate

# Separate header and body
grep "#" ${VCF_FOLDER}/raw.initial2.vcf > ${VCF_FOLDER}/Manta_header
grep -v "#" ${VCF_FOLDER}/raw.initial2.vcf > ${VCF_FOLDER}/Manta_body

# Process MantaINV inversions (INV5/INV3)
cp ${VCF_FOLDER}/Manta_body ${VCF_FOLDER}/Manta_body_edit
cat ${VCF_FOLDER}/Manta_body | grep "MantaINV" | grep "EVENT" | grep "INV5" | cut -f3 | awk -F ':' '{print $1":"$2":"$3":"$4":"}' | while read line ;
do 
    start1=$(grep ${line} ${VCF_FOLDER}/Manta_body | head -1 | cut -f2)
    start2=$(grep ${line} ${VCF_FOLDER}/Manta_body | tail -1 | cut -f2)
    if [[ ${start1} -gt ${start2} ]]
    then
        rmID=$(grep ${line} ${VCF_FOLDER}/Manta_body | tail -1 | cut -f3)
    else
        rmID=$(grep ${line} ${VCF_FOLDER}/Manta_body | head -1 | cut -f3)
    fi
    awk -v rmID=$rmID '$3!=rmID {print $0}' ${VCF_FOLDER}/Manta_body_edit > ${VCF_FOLDER}/tmp 
    mv ${VCF_FOLDER}/tmp ${VCF_FOLDER}/Manta_body_edit
done 

# Process other inversions with 2 partner entries
cp ${VCF_FOLDER}/Manta_body_edit ${VCF_FOLDER}/Manta_body_edit2
grep "INV" ${VCF_FOLDER}/Manta_body_edit | grep -v "EVENT" | cut -f3 | awk -F ':' '{print $1":"$2":"$3":"$4":"}' | while read line ;
do 
    start1=$(grep ${line} ${VCF_FOLDER}/Manta_body_edit | head -1 | cut -f2)
    start2=$(grep ${line} ${VCF_FOLDER}/Manta_body_edit | tail -1 | cut -f2)
    if [[ ${start1} -gt ${start2} ]]
    then
        rmID=$(grep ${line} ${VCF_FOLDER}/Manta_body_edit | tail -1 | cut -f3)
        awk -v rmID=$rmID '$3!=rmID {print $0}' ${VCF_FOLDER}/Manta_body_edit2 > ${VCF_FOLDER}/tmp 
        mv ${VCF_FOLDER}/tmp ${VCF_FOLDER}/Manta_body_edit2
    elif [[ ${start1} -lt ${start2} ]]	
    then
        rmID=$(grep ${line} ${VCF_FOLDER}/Manta_body_edit | head -1 | cut -f3)
        awk -v rmID=$rmID '$3!=rmID {print $0}' ${VCF_FOLDER}/Manta_body_edit2 > ${VCF_FOLDER}/tmp 
        mv ${VCF_FOLDER}/tmp ${VCF_FOLDER}/Manta_body_edit2
    fi
done

# Reconstruct VCF
cat ${VCF_FOLDER}/Manta_header ${VCF_FOLDER}/Manta_body_edit2 > ${VCF_FOLDER}/raw.initial3.vcf

# Sort VCF
module purge
module load BCFtools/1.9-intel-2018b
bcftools sort -Ov -o ${VCF_FOLDER}/raw.initial3.sorted.vcf ${VCF_FOLDER}/raw.initial3.vcf

# Filter out TRA & BND
bcftools filter -i'INFO/SVTYPE!="TRA" & INFO/SVTYPE!="BND"' -o ${VCF_FOLDER}/raw_sorted.noTRA.vcf -Ov ${VCF_FOLDER}/raw.initial3.sorted.vcf
echo "Total number of SVs restricted to INS, DEL, INV, DUP:"
grep -v "#" ${VCF_FOLDER}/raw_sorted.noTRA.vcf | wc -l 

# Filter for SVs >= 50bp
bcftools filter -i'INFO/SVLEN>=50 | INFO/SVLEN<=-50' -o ${VCF_FOLDER}/raw_sorted.noTRA_50bp.vcf -Ov ${VCF_FOLDER}/raw_sorted.noTRA.vcf
echo "Total number of SVs > 50bp:"
grep -v "#" ${VCF_FOLDER}/raw_sorted.noTRA_50bp.vcf | wc -l 

# Add explicit sequences using R script (if available)
# Note: You'll need to adapt the R script path or skip this step if not available
# For now, we'll proceed without it since it's sample-specific
module purge
module load R/4.1.0-foss-2021a
# If you have the R script:
# Rscript /path/to/add_explicit_seq_manta.r ${VCF_FOLDER}/raw_sorted.noTRA_50bp.vcf ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ.vcf ${REFERENCE}

# If R script not available, use the file as-is
cp ${VCF_FOLDER}/raw_sorted.noTRA_50bp.vcf ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ.vcf

module purge
module load HTSlib/1.9-intel-2018b
module load BCFtools/1.9-intel-2018b
module load VCFtools/0.1.16-intel-2018b-Perl-5.28.0

# Export sequences for advanced filtering
bcftools query -f '%CHROM %POS %INFO/END %INFO/SVTYPE %INFO/SVLEN %REF %ALT\n' ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ.vcf > ${VCF_FOLDER}/SV_data_with_seq.txt

# Create blacklists
# Blacklist: N string > 10
grep -P "N{10,}" ${VCF_FOLDER}/SV_data_with_seq.txt | awk '{print $1 "\t" $2 "\t" $6 "\t" $7}' > ${VCF_FOLDER}/N10_blacklist.bed
echo "SVs excluded because of >10N:" 
wc -l ${VCF_FOLDER}/N10_blacklist.bed 

# Blacklist: missing ref sequence
cat ${VCF_FOLDER}/SV_data_with_seq.txt | awk '{if ($6 == "N") print $1 "\t" $2 "\t" $6 "\t" $7;}' > ${VCF_FOLDER}/N_blacklist.bed
echo "SVs excluded because absence of sequence ref:" 
wc -l ${VCF_FOLDER}/N_blacklist.bed 

# Blacklist: missing alt sequence
cat ${VCF_FOLDER}/SV_data_with_seq.txt | awk '{print $1 "\t" $2 "\t" $6 "\t" $7}' | grep -P "<" > ${VCF_FOLDER}/N_blacklist_bis.bed
echo "SVs excluded because absence of sequence alt:" 
wc -l ${VCF_FOLDER}/N_blacklist_bis.bed 

# Combine blacklists
cat ${VCF_FOLDER}/N_blacklist.bed ${VCF_FOLDER}/N_blacklist_bis.bed ${VCF_FOLDER}/N10_blacklist.bed | sort -k1,1 -k2,2n > ${VCF_FOLDER}/blacklist.bed
bgzip -c ${VCF_FOLDER}/blacklist.bed > ${VCF_FOLDER}/blacklist.bed.gz
tabix -s1 -b2 -e2 ${VCF_FOLDER}/blacklist.bed.gz

# Remove blacklisted variants
bcftools view -T ^${VCF_FOLDER}/blacklist.bed.gz ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ.vcf > ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ_Nfiltered.vcf
echo "SVs after filtration for N sequences:" 
grep -v "#" ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ_Nfiltered.vcf | wc -l 

# Final sort
bcftools sort -o ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ_Nfiltered.sorted.vcf -Ov ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ_Nfiltered.vcf

# Keep the filtered VCF
cp ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ_Nfiltered.sorted.vcf ${OUTPUT_VCF}
bgzip -c ${OUTPUT_VCF} > ${OUTPUT_VCF}.gz
tabix ${OUTPUT_VCF}.gz

echo "Filtering complete for ${SAMPLE_NAME}"
echo "Final output: ${OUTPUT_VCF}.gz"
