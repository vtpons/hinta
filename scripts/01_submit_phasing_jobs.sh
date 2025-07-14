#!/bin/bash
#SBATCH --job-name=phasing_job
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=3gb

ROOT_PATH=/path/to/project

for i in {1..22}; do
  echo "#!/bin/bash
#SBATCH --job-name=phasing_chr$i
#SBATCH --output=${ROOT_PATH}/phasing/logs/phasing_chr$i.out
#SBATCH --error=${ROOT_PATH}/phasing/logs/phasing_chr$i.err
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=70gb


${ROOT_PATH}/software/shapeit5/phase_common_static --region $i \\
  --input ${ROOT_PATH}/data/input/chr$i.input.vcf.gz \\
  --pedigree ${ROOT_PATH}/data/input/offspring.fam \\
  --map ${ROOT_PATH}/genetic_maps/b37/chr$i.b37.gmap \\
  --output ${ROOT_PATH}/phasing/output/chr$i.phased.bcf \\
  --thread 8" > run_shapeit5_chr$i.sh

  sbatch run_shapeit5_chr$i.sh
done