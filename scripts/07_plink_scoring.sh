#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb

module load PLINK/1.9-beta6-20190617

INPUT_DIR=/path/to/genomewide/files
WEIGHTS_FILE=/path/to/weights_file.txt
OUTPUT_DIR=/path/to/output

# Score transmitted
plink --bfile ${INPUT_DIR}/genomewide_transmitted \
  --score ${WEIGHTS_FILE} 1 2 3 header sum \
  --out ${OUTPUT_DIR}/transmitted_scores

# Score non-transmitted
plink --bfile ${INPUT_DIR}/genomewide_nontransmitted \
  --score ${WEIGHTS_FILE} 1 2 3 header sum \
  --out ${OUTPUT_DIR}/nontransmitted_scores
