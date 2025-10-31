#!/bin/bash
#SBATCH --job-name=build_Lumpy
#SBATCH --output=/storage/homefs/kw23y068/logfiles/install_Lumpy_%j.log
#SBATCH --error=/storage/homefs/kw23y068/logfiles/install_Lumpy_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=80G
#SBATCH --time=10:00:00
#SBATCH --partition=epyc2

# Exit immediately if a command fails
set -e

# Load compatible modules
#module load GCCcore/10.2.0
module load GCC/10.3.0
#module load CMake/3.18.4-GCCcore-10.2.0
module load CMake/3.20.1-GCCcore-10.3.0

echo "=== Compiler Information ==="
which gcc
gcc --version
which g++
g++ --version

# Compiler setup
export CC=$(which gcc)
export CXX=$(which g++)
export CFLAGS="-O2"
export CXXFLAGS="-O2 -std=c++11"

# Force use of standard ld instead of gold linker
export LDFLAGS="-fuse-ld=bfd"

# Ensure library paths are set
export LD_LIBRARY_PATH=/software.9/software/GCCcore/10.3.0/lib64:$LD_LIBRARY_PATH
export LIBRARY_PATH=/software.9/software/GCCcore/10.3.0/lib64:$LIBRARY_PATH

# Directories
SOFTWARE_DIR="/storage/homefs/kw23y068/software"
LUMPY_DIR="${SOFTWARE_DIR}/Lumpy"
mkdir -p $LUMPY_DIR
cd $LUMPY_DIR

# Clean previous attempts
rm -rf ./lumpy-sv

# Clone LUMPY repository including submodules
git clone --recursive https://github.com/arq5x/lumpy-sv.git
cd lumpy-sv

# Build LUMPY (verbose)
make clean || true
make -j2 V=1

# Install binaries locally
mkdir -p $LUMPY_DIR/bin
cp bin/* $LUMPY_DIR/bin

# Add LUMPY to PATH temporarily for testing
export PATH=$LUMPY_DIR/bin:$PATH
echo "=== LUMPY installation complete ==="
lumpyexpress --version || lumpy --version

# Instructions for permanent usage
echo "LUMPY installed to $LUMPY_DIR/bin"
echo "To use it permanently, add the following line to your ~/.bashrc:"
echo "export PATH=$LUMPY_DIR/bin:\$PATH"
