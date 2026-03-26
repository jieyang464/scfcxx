#!/bin/bash
# Downloads basis set files for libint2 from the official repository
# and places them where libint2 expects them (/usr/share/libint/2.7.1/basis/)
# or alternatively under LIBINT_DATA_PATH

set -e

# Choose where to install basis sets
# Option A: System location (requires sudo)
# BASIS_DIR="/usr/share/libint/2.7.1/basis"
# sudo mkdir -p "$BASIS_DIR"

# Option B: Under your project (no sudo needed)
BASIS_DIR="$(dirname "$0")/third_party/libint2_basis"
mkdir -p "$BASIS_DIR"

echo "Installing basis sets to: $BASIS_DIR"

# Download from the official libint repository (basis files are under export/lib/basis)
BASE_URL="https://raw.githubusercontent.com/evaleev/libint/master/export/lib/basis"

# Common basis sets
BASIS_SETS=(
    "sto-3g"
    "3-21g"
    "6-31g"
    "6-31g_st_"        # 6-31G*
    "6-31g_st__st_"    # 6-31G**
    "6-31_pp_g"        # 6-31+G
    "6-311g"
    "cc-pvdz"
    "cc-pvtz"
    "cc-pvqz"
    "aug-cc-pvdz"
    "aug-cc-pvtz"
    "def2-svp"
    "def2-tzvp"
    "def2-tzvpp"
    "def2-qzvpp"
)

for basis in "${BASIS_SETS[@]}"; do
    echo "Downloading ${basis}.g94 ..."
    curl -sL "${BASE_URL}/${basis}.g94" -o "${BASIS_DIR}/${basis}.g94" 2>/dev/null || \
    wget -q "${BASE_URL}/${basis}.g94" -O "${BASIS_DIR}/${basis}.g94" 2>/dev/null || \
    echo "  WARNING: Failed to download ${basis}.g94"
done

echo ""
echo "Done! Basis sets installed to: $BASIS_DIR"
echo ""
echo "To use them, set the environment variable before running your program:"
echo "  export LIBINT_DATA_PATH=${BASIS_DIR}"
echo ""
echo "Or add this to your .bashrc:"
echo "  echo 'export LIBINT_DATA_PATH=${BASIS_DIR}' >> ~/.bashrc"
