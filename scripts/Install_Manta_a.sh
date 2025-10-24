#!/bin/bash
#SBATCH --job-name=install_manta_prebuilt
#SBATCH --output=/storage/homefs/kw23y068/logfiles/install_manta_prebuilt_%j.out
#SBATCH --error=/storage/homefs/kw23y068/logfiles/install_manta_prebuilt_%j.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

echo "=== Installing pre-built Manta ==="
echo "Date: $(date)"

cd /storage/homefs/kw23y068/software/

# Remove old failed installation
echo ">>> Cleaning up old installation"
rm -rf Manta

# Download pre-built binary
echo ">>> Downloading pre-built Manta"
wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2

# Extract
echo ">>> Extracting"
tar -xjf manta-1.6.0.centos6_x86_64.tar.bz2

# Create symlink for easy access
ln -sf manta-1.6.0.centos6_x86_64 Manta

# Test it works
echo ">>> Testing installation"
./Manta/bin/configManta.py --version

echo "=== Installation complete ==="
echo "Manta installed at: /storage/homefs/kw23y068/software/Manta"
echo "Date: $(date)"