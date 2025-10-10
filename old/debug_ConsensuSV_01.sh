#!/bin/bash
#SBATCH --job-name=consensusv_debug
#SBATCH --output=../logfiles/consensusv_debug_%j.out
#SBATCH --error=../logfiles/consensusv_debug_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=epyc2

cd /storage/homefs/kw23y068/software/ConsensuSV-pipeline


# Run the debug command
apptainer exec \
  --bind /storage/homefs/kw23y068/software/ConsensuSV-pipeline:/ConsensuSV-pipeline \
  --bind /storage/homefs/kw23y068/test_working_dir:/test_working_dir \
  --bind /storage/homefs/kw23y068/tools:/tools \
  --bind /storage/homefs/kwe23y068/logfiles/luigi:/logfiles/luigi \
  /storage/homefs/kw23y068/software/consensusv-pipeline.sif \
  bash -c "python -u /tools/ConsensuSV-core/main.py \
  -of /test_working_dir/output/ \
  -f /test_working_dir/pipeline/ \
  -s ERR018471 \
  -c breakdancer,breakseq,cnvnator,delly,lumpy,manta,tardis,whamg \
  -mod /tools/ConsensuSV-core/pretrained_1000g_illumina.model"
