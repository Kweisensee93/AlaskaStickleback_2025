#!/bin/bash
#SBATCH --job-name=install_manta
#SBATCH --output=/storage/homefs/kw23y068/logfiles/install_manta_%j.log
#SBATCH --error=/storage/homefs/kw23y068/logfiles/install_manta_%j.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G

echo "=== Starting Manta installation job on $(hostname) ==="
echo "Date: $(date)"

### --- Environment check ---
echo ">>> Checking compiler and Python versions"
gcc --version
python --version
module load CMake/3.29.3-GCCcore-13.3.0

### --- Set up directories ---
export BASEDIR=/storage/homefs/kw23y068/software/Manta
export SRC=$BASEDIR/src
export BUILD=$BASEDIR/build
export INSTALL=$BASEDIR/install
# Check before installation, if newer is available
export MANTA_VERSION=1.6.0


if [ ! -d "$SRC" ]
then
    mkdir $SRC
fi
if [ ! -d "$BUILD" ]
then
    mkdir $BUILD
fi
if [ ! -d "$INSTALL" ]
then
    mkdir $INSTALL
fi

cd $SRC

### --- Download source if missing ---
if [ ! -f manta-${MANTA_VERSION}.release_src.tar.bz2 ]; then
    echo ">>> Downloading Manta v${MANTA_VERSION}"
    wget https://github.com/Illumina/manta/releases/download/v${MANTA_VERSION}/manta-${MANTA_VERSION}.release_src.tar.bz2
fi

### --- Unpack source ---
tar -xjf manta-${MANTA_VERSION}.release_src.tar.bz2

### --- Build Manta ---
cd $BUILD

../src/manta-${MANTA_VERSION}.release_src/configure \
    --jobs=${SLURM_CPUS_PER_TASK:-4} \
    --prefix=$INSTALL

echo ">>> Installation complete"
echo "Installed at: $INSTALL"

echo "=== Manta installation job completed ==="
echo "Date: $(date)"
