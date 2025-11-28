#!/bin/bash
#SBATCH --job-name=lumpy_filter
#SBATCH --output=/storage/homefs/kw23y068/logfiles/lumpy_filter_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/lumpy_filter_%j.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=epyc2

## This script is based on a smoove script.
## smoove runs lumpy but also adds some filtering

####################
# PARAMETERS TO SET
####################
# May be altered for feeding in arguments to script
# Define which batch of samples to process
BATCH_START=1    # Change to 81, 161, 241, 321, 401 for other batches
BATCH_END=80     # Change to 160, 240, 320, 400, 480 for other batches

START_PAD=$(printf "%03d" "${BATCH_START}")
END_PAD=$(printf "%03d" "${BATCH_END}")
INPUT_DIR="/storage/scratch/iee_evol/kw23y068/lumpy_joint_${START_PAD}_${END_PAD}"



gunzip SPI_BJO_ROS-smoove.genotyped.vcf.gz 
cp SPI_BJO_ROS-smoove.genotyped.vcf smoove_merged.vcf

echo "nb of SV detected by Smoove"
grep -v "#" smoove_merged.vcf | wc -l 

module purge
module load HTSlib/1.9-intel-2018b
module load BCFtools/1.9-intel-2018b
module load VCFtools/0.1.16-intel-2018b-Perl-5.28.0

bgzip smoove_merged.vcf

#Filter by developers
bcftools view -i '(SVTYPE = "DEL" & FMT/DHFFC[0] < 0.7) | (SVTYPE = "DUP" & FMT/DHBFC[0] > 1.3) | (SVTYPE = "INV") | (SVTYPE = "INS")' -O z -o smoove_merged.filter1.vcf.gz smoove_merged.vcf.gz
bcftools sort -O z -o smoove_merged.filter1.sorted.vcf.gz smoove_merged.filter1.vcf.gz ;

zcat smoove_merged.filter1.sorted.vcf.gz | grep -v "#" | wc -l 

INPUT_VCF=smoove_merged.filter1.sorted.vcf.gz  # we will try not to modify this one
VCF_FOLDER=TMP_VCF
OUTPUT_VCF=SPI_BJO_ROS_smoove.vcf
RScripts=SV_Merot/Rscripts
RefFasta=SV_Merot/Ref/Puffin_ReferenceGenome.NU.MT.fasta

mkdir ${VCF_FOLDER}
cp $INPUT_VCF ${VCF_FOLDER}/raw.vcf.gz
gunzip ${VCF_FOLDER}/raw.vcf.gz

#filter vcf -i (include, -O vcf format -o
bcftools filter -i'INFO/SVLEN>=50 | INFO/SVLEN<=-50' -o ${VCF_FOLDER}/raw_sorted.noTRA_50bp.vcf -Ov ${VCF_FOLDER}/raw.vcf
echo "total number of SVs > 50b"
grep -v "#" ${VCF_FOLDER}/raw_sorted.noTRA_50bp.vcf | wc -l #2,141

#then we use the reference to get the sequence with a R scripts graciously provided by Marc-André Lemay
module purge
module load R/4.1.0-foss-2021a
Rscript ${RScripts}/add_explicit_seq.r ${VCF_FOLDER}/raw_sorted.noTRA_50bp.vcf ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ.vcf ${RefFasta}

module purge
module load HTSlib/1.9-intel-2018b
module load BCFtools/1.9-intel-2018b
module load VCFtools/0.1.16-intel-2018b-Perl-5.28.0
#Export sequences for advanced filtering
bcftools query -f '%CHROM %POS %INFO/END %INFO/SVTYPE %INFO/SVLEN %REF %ALT\n' ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ.vcf > ${VCF_FOLDER}/SV_data_with_seq.txt

#blacklist because of N string > 10 (possible junction of contigs 
grep -P "N{10,}" ${VCF_FOLDER}/SV_data_with_seq.txt | awk '{print $1 "\t" $2 "\t" $6 "\t" $7}' > ${VCF_FOLDER}/N10_blacklist.bed
echo "SVs excluded because of >10N" 
wc -l ${VCF_FOLDER}/N10_blacklist.bed 

#blacklist because missing seq
cat ${VCF_FOLDER}/SV_data_with_seq.txt | awk '{if ($6 == "N") print $1 "\t" $2 "\t" $6 "\t" $7;}' > ${VCF_FOLDER}/N_blacklist.bed
echo "SVs excluded because absence of sequence ref" 
wc -l ${VCF_FOLDER}/N_blacklist.bed 

#blacklist because missing seq
cat  ${VCF_FOLDER}/SV_data_with_seq.txt | awk '{print $1 "\t" $2 "\t" $6 "\t" $7}' | grep -P "<" > ${VCF_FOLDER}/N_blacklist_bis.bed
echo "SVs excluded because absence of sequence alt" 
wc -l ${VCF_FOLDER}/N_blacklist_bis.bed 

#full blacklist
cat ${VCF_FOLDER}/N_blacklist.bed ${VCF_FOLDER}/N_blacklist_bis.bed ${VCF_FOLDER}/N10_blacklist.bed | sort -k1,1 -k2,2n > ${VCF_FOLDER}/blacklist.bed
#head $VCF_FOLDER/blacklist.bed
bgzip -c ${VCF_FOLDER}/blacklist.bed > ${VCF_FOLDER}/blacklist.bed.gz
tabix -s1 -b2 -e2 ${VCF_FOLDER}/blacklist.bed.gz

#remove blacklist of variants
bcftools view -T ^${VCF_FOLDER}/blacklist.bed.gz ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ.vcf > ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ_Nfiltered.vcf
echo "SVs after filtration for N seq" 
grep -v "#" ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ_Nfiltered.vcf | wc -l 

#keep the filtered vcf
cp ${VCF_FOLDER}/raw_sorted.noTRA_50bp_withSEQ_Nfiltered.vcf ${OUTPUT_VCF}
bgzip -c ${OUTPUT_VCF} > ${OUTPUT_VCF}.gz
tabix ${OUTPUT_VCF}.gz
