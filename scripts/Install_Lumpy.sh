#!/bin/bash
#SBATCH --job-name=build_Lumpy
#SBATCH --output=/storage/homefs/kw23y068/logfiles/install_Lumpy_%j.log
#SBATCH --error=/storage/homefs/kw23y068/logfiles/install_Lumpy_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=80G
#SBATCH --time=10:00:00
#SBATCH --partition=epyc2

# Exit on any error
set -e

# Prerequesites according to Github:
# G++ compiler (is pre-installed on the cluster)
# CMake
# module load CMake/3.29.3-GCCcore-13.3.0
## due to compatibility, we use older CMake
## module load CMake/3.21.1-GCCcore-11.2.0
# Incompatibility is due to too new C++ not CMake
# export CXXFLAGS="-std=c++11"
# export CFLAGS="-std=c11"
module load GCCcore/10.2.0
module load CMake/3.18.4-GCCcore-10.2.0

export CC=$(which gcc)
export CXX=$(which g++)


# Directories
SOFTWARE_DIR="/storage/homefs/kw23y068/software"
mkdir -p ${SOFTWARE_DIR}/Lumpy
cd ${SOFTWARE_DIR}/Lumpy
# clean previous attempt
rm -rf ./lumpy-sv

# Clone and build LUMPY
git clone --recursive https://github.com/arq5x/lumpy-sv.git
cd lumpy-sv
make clean
make

# Install binaries locally (no admin rights, which would be necessary for GitHub instructions)
mkdir -p ${SOFTWARE_DIR}/Lumpy/bin
cp bin/* ${SOFTWARE_DIR}/Lumpy/bin

# Add to PATH temporarily for testing
export PATH=${LUMPY_DIR}/bin:$PATH
echo "=== LUMPY installation complete ==="
lumpyexpress --version || lumpy --version

echo "LUMPY installed to ${SOFTWARE_DIR}/Lumpy/bin"
echo "To use it, add the following line to your ~/.bashrc:"
echo "export PATH=${SOFTWARE_DIR}/Lumpy/bin:\$PATH"
echo "or manually add it to scripts as ${SOFTWARE_DIR}/Lumpy/bin"
