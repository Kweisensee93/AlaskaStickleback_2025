#!/bin/bash
#SBATCH --job-name=consensusv_test
#SBATCH --output=../logfiles/consensusv_test_%j.out
#SBATCH --error=../logfiles/consensusv_test_%j.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=epyc2

# Ensure working dir exists
mkdir -p /storage/homefs/kw23y068/test_working_dir
mkdir -p /storage/homefs/kw23y068/logfiles/luigi

export TMPDIR=/storage/homefs/kw23y068/test_working_dir/tmp
mkdir -p $TMPDIR


cd /storage/homefs/kw23y068/software/ConsensuSV-pipeline

# Start luigid with a writable logdir
#apptainer exec /storage/homefs/kw23y068/software/consensusv-pipeline.sif \
#    luigid --background --port 8082 --logdir /storage/homefs/kw23y068/logfiles/luigi

# Run the test
# We need a test woring directory with write permissions because we use apptainer insted of docker
# we bind whole /storage/homefs/kw23y068 to /home to ensure that all paths in the test script are correct
# we bind /storage/homefs/kw23y068/logfiles/luigi for logfiles
# we bind /storage/homefs/kw23y068/test_working_dir for tmp files
# luigi-id and test rund need to be started in the same apptainer instance
# Everything is separated by && to ensure that the next command is only executed if the previous one was successful
# The sleep 5s may be redundant to the && commands, but it ensures that luigid has enough time to start before the test run is executed
apptainer exec \
    --bind /storage/homefs/kw23y068:/home \
    --bind /storage/homefs/kw23y068/test_working_dir:/test_working_dir \
    --bind /storage/homefs/kw23y068/logfiles/luigi:/logfiles/luigi \
    /storage/homefs/kw23y068/software/consensusv-pipeline.sif \
    bash -c \
    "luigid --background --port 8082 --logdir /home/logfiles/luigi \
    && sleep 5s \
    && cd /home/software/ConsensuSV-pipeline \
    && /home/software/ConsensuSV-pipeline/test_run_csv.sh"


# possible optimization for the sleep time:
# until nc -z localhost 8082; do
#         sleep 1

# as for now we skip multiple bindings and simplify to /storage/homefs/kw23y068:/home
# binds kicked out:
# --bind /storage/homefs/kw23y068/software/ConsensuSV-pipeline:/ConsensuSV-pipeline \
# --bind /storage/homefs/kw23y068/test_working_dir:/test_working_dir \
# --bind /storage/homefs/kw23y068/logfiles/luigi:/logfiles/luigi \
# --bind /storage/homefs/kw23y068/test_working_dir:/test_working_dir \

# we have a fail in finding stuff:
    #--bind /storage/homefs/kw23y068/software/ConsensuSV-core:/tools/ConsensuSV-core \
# this does not fix any of this