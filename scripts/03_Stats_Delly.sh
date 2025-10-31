#!/bin/bash


module load BCFtools/1.12-GCC-10.3.0


# Script to generate status.txt files for already-processed Delly samples
# This checks the output files to determine if the run was successful

DELLY_BASE_DIR="/storage/scratch/iee_evol/kw23y068/Delly"

echo "=== Generating status files for existing Delly results ==="
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
        # Don't skip - will overwrite
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
        # Still count as failed since pipeline didn't complete
        ((FAILED_COUNT++))
    else
        # All files exist - check if they have content
        BCF_SIZE=$(stat -f%z "${BCF_FILE}" 2>/dev/null || stat -c%s "${BCF_FILE}" 2>/dev/null)
        
        if [ "${BCF_SIZE}" -lt 100 ]; then
            echo "FAILED: BCF too small (${BCF_SIZE} bytes)" > "${STATUS_FILE}"
            echo "  -> BCF suspiciously small"
            ((FAILED_COUNT++))
        else
            # Success! Count SVs if bcftools is available
            if command -v bcftools &> /dev/null; then
                TOTAL_SVS=$(bcftools view -H "${VCF_GZ}" 2>/dev/null | wc -l)
                PASS_SVS=$(bcftools view -H "${FILTERED_VCF}" 2>/dev/null | wc -l)
                
                echo "SUCCESS" > "${STATUS_FILE}"
                echo "Total SVs: ${TOTAL_SVS}" >> "${STATUS_FILE}"
                echo "PASS SVs: ${PASS_SVS}" >> "${STATUS_FILE}"
                echo "  -> Success (Total: ${TOTAL_SVS}, PASS: ${PASS_SVS})"
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
echo "Skipped (already had status): see above"
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
