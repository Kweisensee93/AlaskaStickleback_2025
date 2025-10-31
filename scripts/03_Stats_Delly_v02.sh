#!/bin/bash
#SBATCH --job-name=delly_stats
#SBATCH --output=/storage/homefs/kw23y068/logfiles/delly_stats_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/delly_stats_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --partition=epyc2

# Script to generate detailed status.txt files for already-processed Delly samples
# Includes SV type breakdowns for both original and filtered VCF files

module load BCFtools/1.12-GCC-10.3.0

DELLY_BASE_DIR="/storage/scratch/iee_evol/kw23y068/Delly"

echo "=== Generating detailed status files for existing Delly results ==="
echo "Date: $(date)"
echo ""

# Counters
SUCCESS_COUNT=0
FAILED_COUNT=0
SKIPPED_COUNT=0
TOTAL_COUNT=0

# Loop through all sample directories
for SAMPLE_DIR in ${DELLY_BASE_DIR}/*/; do
    if [ ! -d "${SAMPLE_DIR}" ]; then
        continue
    fi
    
    SAMPLE_NAME=$(basename "${SAMPLE_DIR}")
    STATUS_FILE="${SAMPLE_DIR}/status.txt"
    
    if [ -f "${STATUS_FILE}" ]; then
        echo "Overwriting existing status file for ${SAMPLE_NAME}..."
    fi
    
    ((TOTAL_COUNT++))
    
    # Define expected files
    BCF_FILE="${SAMPLE_DIR}/${SAMPLE_NAME}.bcf"
    VCF_FILE="${SAMPLE_DIR}/${SAMPLE_NAME}.vcf"
    VCF_GZ="${SAMPLE_DIR}/${SAMPLE_NAME}.vcf.gz"
    FILTERED_VCF="${SAMPLE_DIR}/${SAMPLE_NAME}.filtered.vcf.gz"
    
    echo "Checking ${SAMPLE_NAME}..."
    
    # Check for various failure scenarios
    if [ ! -f "${BCF_FILE}" ]; then
        echo "FAILED: No BCF output" > "${STATUS_FILE}"
        echo "  -> No BCF file found"
        ((FAILED_COUNT++))
    elif [ ! -s "${BCF_FILE}" ]; then
        echo "FAILED: Empty BCF" > "${STATUS_FILE}"
        echo "  -> BCF file is empty"
        ((FAILED_COUNT++))
    elif [ ! -f "${VCF_GZ}" ]; then
        echo "FAILED: No VCF output" > "${STATUS_FILE}"
        echo "  -> No compressed VCF found"
        ((FAILED_COUNT++))
    elif [ ! -f "${FILTERED_VCF}" ]; then
        echo "PARTIAL: Missing filtered VCF" > "${STATUS_FILE}"
        echo "  -> Filtered VCF missing, but main files exist"
        ((FAILED_COUNT++))
    else
        # All files exist - check if they have content
        BCF_SIZE=$(stat -f%z "${BCF_FILE}" 2>/dev/null || stat -c%s "${BCF_FILE}" 2>/dev/null)
        
        if [ "${BCF_SIZE}" -lt 100 ]; then
            echo "FAILED: BCF too small (${BCF_SIZE} bytes)" > "${STATUS_FILE}"
            echo "  -> BCF suspiciously small"
            ((FAILED_COUNT++))
        else
            # Success! Count SVs and get detailed type information
            if command -v bcftools &> /dev/null; then
                # Total SVs in original VCF
                TOTAL_SVS=$(bcftools view -H "${VCF_GZ}" 2>/dev/null | wc -l)
                
                # Total SVs in filtered VCF
                PASS_SVS=$(bcftools view -H "${FILTERED_VCF}" 2>/dev/null | wc -l)
                
                # Count SV types in original VCF
                DEL_COUNT=$(bcftools query -f '%SVTYPE\n' "${VCF_GZ}" 2>/dev/null | grep -c "^DEL$" || echo "0")
                DUP_COUNT=$(bcftools query -f '%SVTYPE\n' "${VCF_GZ}" 2>/dev/null | grep -c "^DUP$" || echo "0")
                INV_COUNT=$(bcftools query -f '%SVTYPE\n' "${VCF_GZ}" 2>/dev/null | grep -c "^INV$" || echo "0")
                BND_COUNT=$(bcftools query -f '%SVTYPE\n' "${VCF_GZ}" 2>/dev/null | grep -c "^BND$" || echo "0")
                INS_COUNT=$(bcftools query -f '%SVTYPE\n' "${VCF_GZ}" 2>/dev/null | grep -c "^INS$" || echo "0")
                
                # Count SV types in filtered VCF
                DEL_PASS=$(bcftools query -f '%SVTYPE\n' "${FILTERED_VCF}" 2>/dev/null | grep -c "^DEL$" || echo "0")
                DUP_PASS=$(bcftools query -f '%SVTYPE\n' "${FILTERED_VCF}" 2>/dev/null | grep -c "^DUP$" || echo "0")
                INV_PASS=$(bcftools query -f '%SVTYPE\n' "${FILTERED_VCF}" 2>/dev/null | grep -c "^INV$" || echo "0")
                BND_PASS=$(bcftools query -f '%SVTYPE\n' "${FILTERED_VCF}" 2>/dev/null | grep -c "^BND$" || echo "0")
                INS_PASS=$(bcftools query -f '%SVTYPE\n' "${FILTERED_VCF}" 2>/dev/null | grep -c "^INS$" || echo "0")
                
                # Write detailed status file
                echo "SUCCESS" > "${STATUS_FILE}"
                echo "" >> "${STATUS_FILE}"
                echo "=== Original VCF (${SAMPLE_NAME}.vcf.gz) ===" >> "${STATUS_FILE}"
                echo "Total SVs: ${TOTAL_SVS}" >> "${STATUS_FILE}"
                echo "  DEL: ${DEL_COUNT}" >> "${STATUS_FILE}"
                echo "  DUP: ${DUP_COUNT}" >> "${STATUS_FILE}"
                echo "  INV: ${INV_COUNT}" >> "${STATUS_FILE}"
                echo "  INS: ${INS_COUNT}" >> "${STATUS_FILE}"
                echo "  BND: ${BND_COUNT}" >> "${STATUS_FILE}"
                echo "" >> "${STATUS_FILE}"
                echo "=== Filtered VCF (${SAMPLE_NAME}.filtered.vcf.gz) ===" >> "${STATUS_FILE}"
                echo "PASS SVs: ${PASS_SVS}" >> "${STATUS_FILE}"
                echo "  DEL: ${DEL_PASS}" >> "${STATUS_FILE}"
                echo "  DUP: ${DUP_PASS}" >> "${STATUS_FILE}"
                echo "  INV: ${INV_PASS}" >> "${STATUS_FILE}"
                echo "  INS: ${INS_PASS}" >> "${STATUS_FILE}"
                echo "  BND: ${BND_PASS}" >> "${STATUS_FILE}"
                
                echo "  -> Success (Total: ${TOTAL_SVS} [DEL:${DEL_COUNT} DUP:${DUP_COUNT} INV:${INV_COUNT} INS:${INS_COUNT} BND:${BND_COUNT}], PASS: ${PASS_SVS} [DEL:${DEL_PASS} DUP:${DUP_PASS} INV:${INV_PASS} INS:${INS_PASS} BND:${BND_PASS}])"
                ((SUCCESS_COUNT++))
            else
                echo "SUCCESS" > "${STATUS_FILE}"
                echo "Total SVs: Unknown (bcftools not loaded)" >> "${STATUS_FILE}"
                echo "PASS SVs: Unknown (bcftools not loaded)" >> "${STATUS_FILE}"
                echo "  -> Success (counts unavailable)"
                ((SUCCESS_COUNT++))
            fi
        fi
    fi
done

echo ""
echo "=== Summary ==="
echo "Total samples processed: ${TOTAL_COUNT}"
echo "Successful: ${SUCCESS_COUNT}"
echo "Failed: ${FAILED_COUNT}"
echo ""
echo "Date: $(date)"
echo ""
echo "To view all statuses:"
echo "  cat ${DELLY_BASE_DIR}/*/status.txt"
echo ""
echo "To find failed samples:"
echo "  grep -l 'FAILED' ${DELLY_BASE_DIR}/*/status.txt"
echo ""
echo "To get summary counts:"
echo "  grep -h '^SUCCESS\\|^FAILED\\|^PARTIAL' ${DELLY_BASE_DIR}/*/status.txt | sort | uniq -c"
echo ""
echo "To view a specific sample's detailed stats:"
echo "  cat ${DELLY_BASE_DIR}/<sample_name>/status.txt"