#!/bin/bash
#SBATCH --job-name=build_manta
#SBATCH --output=/storage/homefs/kw23y068/logfiles/build_manta_%j.log
#SBATCH --error=/storage/homefs/kw23y068/logfiles/build_manta_%j.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G

echo "=== Building Manta on $(hostname) ==="
echo "Date: $(date)"

# Load Python 2.7 FIRST
module purge
module load Python/2.7.18-GCCcore-11.3.0-bare
module load CMake/3.29.3-GCCcore-13.3.0

# Verify Python version
echo ">>> Using Python version:"
python --version
which python

# Clean previous failed build
cd /storage/homefs/kw23y068/software/Manta
rm -rf build install

# Reconfigure with Python 2.7
mkdir build install
cd build
../src/manta-1.6.0.release_src/configure \
    --jobs=${SLURM_CPUS_PER_TASK} \
    --prefix=/storage/homefs/kw23y068/software/Manta/install

# Build and install
echo ">>> Running make install"
make -j ${SLURM_CPUS_PER_TASK} install

echo "=== Build completed ==="
echo "Date: $(date)"