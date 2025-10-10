#!/bin/bash
#SBATCH --job-name=svmeca_test
#SBATCH --cpus-per-task=16
#SBATCH --mem=70G
#SBATCH --time=04:00:00
#SBATCH --output=/storage/homefs/kw23y068/logfiles/%j_svmeca_test1.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/%j_svmeca_test1.err
#SBATCH --partition=epyc2

set -Eeuo pipefail
set -o verbose

# ressources
CPU=$SLURM_CPUS_PER_TASK
RAM=64G     # use slightly less than allocated memory to avoid OOM kills
REFERENCE="GRCh38_full_analysis_set_plus_decoy_hla.fa"

# get the container (workaround for apptainer)
CONTAINER="docker://wembasop/sv-meca:latest"

# create a local writable temp directory for Java/Nextflow
mkdir -p /storage/homefs/kw23y068/software/SV-MeCa_data/test/tmp
mkdir -p /storage/homefs/kw23y068/output/svmeca/test/tmp

cd /storage/homefs/kw23y068/software/SV-MeCa_data/test/input

# get the cram 60min
wget -q ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR398/ERR3988781/HG00514.final.cram	

# get the reference 10min
wget -q ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/$REFERENCE

cd ..

# cram to bam 7min
#docker run -v $(pwd)/input:/input wembasop/sv-meca:latest "/opt/conda/envs/samtools/bin/samtools view -@ $CPU -b -o /input/HG00514.bam -T /input/$REFERENCE /input/HG00514.final.cram"
apptainer exec \
  -B $(pwd)/input:/input \
  $CONTAINER \
  /opt/conda/envs/samtools/bin/samtools view -@ $CPU -b -o /input/HG00514.bam -T /input/$REFERENCE /input/HG00514.final.cram

rm input/HG00514.final.cram

# downsample 16 min
#docker run -v $(pwd)/input:/input wembasop/sv-meca:latest "/opt/conda/envs/align/bin/gatk --java-options "-Xmx$RAM" DownsampleSam -I /input/HG00514.bam -O /input/HG00514.down.bam -P 0.30" 
apptainer exec \
  -B $(pwd)/input:/input \
  $CONTAINER \
  bash -c "/opt/conda/envs/align/bin/gatk --java-options \"-Xmx$RAM\" DownsampleSam -I /input/HG00514.bam -O /input/HG00514.down.bam -P 0.30"

rm input/HG00514.bam

# test if /tmp is writable
#apptainer exec $CONTAINER bash -c "touch /tmp/testfile && echo tmp_is_writable"
## The .out file shows tmp_is_writable

## For fixing, after 31867090.err
## It needs Java for nextflow and a writable /tmp directory
## Java needs to find its way inside the last apptainer

# run SV-MeCa (141min)
#docker run -v $(pwd)/input:/input -v $(pwd)/output:/workspace/SV-MeCa/results wembasop/sv-meca:latest "/workspace/SV-MeCa/run_svmeca.sh bam -bam /input/HG00514.down.bam -ref /input/$REFERENCE -sample HG00514 -build hg38 -has_chr true -bed /input/hg38_centromer.bed" 
# run SV-MeCa (~141 min) with proper temp directory
apptainer exec \
  -B /storage/homefs/kw23y068/software/SV-MeCa_data/test/input:/input \
  -B /storage/homefs/kw23y068/output/svmeca/test:/workspace/SV-MeCa/results \
  -B /storage/homefs/kw23y068/output/svmeca/test/tmp:/tmp \
  $CONTAINER \
  bash -c "export TMPDIR=/tmp; export NXF_TMP=/tmp; export JAVA_CMD=/opt/conda/envs/align/bin/java; \
  /workspace/SV-MeCa/run_svmeca.sh \
  bam -bam /input/HG00514.down.bam \
  -ref /input/$REFERENCE -sample HG00514 -build hg38 -has_chr true -bed /input/hg38_centromer.bed"

# with #bash -c echo "Hello world, I made it this far!"\ a correct setup was confirmed.
# The export commands for TMDIR and NXF_TMP are needed to make the image work with apptainer
# Within the container CAPSULE fails to autodetect Java, so it is set to /opt/conda/envs/align/bin/java


# clean up temp directory
rm -rf /storage/homefs/kw23y068/output/svmeca/test/tmp/*
