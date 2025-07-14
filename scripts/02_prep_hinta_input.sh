#!/bin/bash
#SBATCH --job-name=prep_hinta_input
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=3gb

ROOT_PATH=/path/to/project

module load PLINK/2.0
module load BCFtools
module load PythonPlus

# Paths 
PHASED_DIR=${ROOT_PATH}/phasing/output
HINTA_INPUT_DIR=${ROOT_PATH}/hinta/input
SCRIPTS_DIR=${ROOT_PATH}/scripts/hinta
MAP_FILE_DIR=${ROOT_PATH}/genetic_maps/b37

# Sample file: can use the same for all iomosomes 
SAMPLE_FILE=${HINTA_INPUT_DIR}/chr22.sample

for i in {1..22}; do

  # Convert BCF to VCF
  bcftools convert \
    -Oz \
    -o ${PHASED_DIR}/i${i}.phased.vcf.gz \
    ${PHASED_DIR}/i${i}.phased.bcf

  bcftools index -f ${PHASED_DIR}/i${i}.phased.vcf.gz

  # Convert VCF to haps/sample
  plink2 \
    --vcf ${PHASED_DIR}/i${i}.phased.vcf.gz \
    --export haps \
    --out ${HINTA_INPUT_DIR}/i${i}.phased

  # Create subsets (sample files with N=1000 each) 
  python ${SCRIPTS_DIR}/subset_phased_data.py \
    ${i} \
    ${HINTA_INPUT_DIR}/i${i}.phased.haps \
    ${SAMPLE_FILE} \
    ${HINTA_INPUT_DIR}/subsets/
done
