#!/bin/bash
#SBATCH --job-name=create_genomewide_PGS_inputs
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=50gb
#SBATCH --nodes=1

module load PLINK/2.0-alpha4.5-20230813
module load BCFtools

ROOT_PATH=/path/to/working_directory

TRANSMITTED_INPUT=${ROOT_PATH}/transmitted/input
TRANSMITTED_OUTPUT=${ROOT_PATH}/transmitted/output
TRANSMITTED_MERGELIST=${TRANSMITTED_INPUT}/mergelist_transmitted.txt
TRANSMITTED_VCF=${TRANSMITTED_INPUT}/genomewide_transmitted.vcf.gz
TRANSMITTED_PREFIX=${TRANSMITTED_OUTPUT}/genomewide_transmitted

NONTRANSMITTED_INPUT=${ROOT_PATH}/nontransmitted/input
NONTRANSMITTED_OUTPUT=${ROOT_PATH}/nontransmitted/output
NONTRANSMITTED_MERGELIST=${NONTRANSMITTED_INPUT}/mergelist_nontransmitted.txt
NONTRANSMITTED_VCF=${NONTRANSMITTED_INPUT}/genomewide_nontransmitted.vcf.gz
NONTRANSMITTED_PREFIX=${NONTRANSMITTED_OUTPUT}/genomewide_nontransmitted

# ======= TRANSMITTED ======= #

# Create merge list 
for i in {1..22}; do 
     echo "${TRANSMITTED_INPUT}/chr${i}.transmitted.split.vcf.gz"
 done > "${TRANSMITTED_MERGELIST}"

# Concatenate and index (optional)
bcftools concat --file-list "${TRANSMITTED_MERGELIST}" --output-type z --output "${TRANSMITTED_VCF}"
bcftools index "${TRANSMITTED_VCF}"

# Convert to PLINK format
plink2 --vcf "${TRANSMITTED_VCF}" --make-bed --double-id --out "${TRANSMITTED_PREFIX}"

# ======= NON-TRANSMITTED ======= # 

# Create merge list
for i in {1..22}; do 
    echo "${NONTRANSMITTED_INPUT}/chr${i}.nontransmitted.split.vcf.gz"
done > "${NONTRANSMITTED_MERGELIST}"

# Concatenate and index 
bcftools concat --file-list "${NONTRANSMITTED_MERGELIST}" --output-type z --output "${NONTRANSMITTED_VCF}"
bcftools index "${NONTRANSMITTED_VCF}"

# Convert to PLINK format
plink2 --vcf "${NONTRANSMITTED_VCF}" --make-bed --double-id --out "${NONTRANSMITTED_PREFIX}"
