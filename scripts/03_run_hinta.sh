#!/bin/bash
#SBATCH --job-name=run_hinta
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb
#SBATCH --nodes=1

ROOT_PATH=/path/to/project

for i in {1..22}
do
  echo "#!/bin/bash
#SBATCH --job-name=run_hinta_chr$i
#SBATCH --output=${ROOT_PATH}/logs/run_hinta_chr$i.out
#SBATCH --error=${ROOT_PATH}/logs/run_hinta_chr$i.err
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb
#SBATCH --nodes=1

ROOT_PATH=/path/to/project

for c in {0..64} 
do
  echo \"#!/bin/bash
#SBATCH --job-name=run_hinta_chr$i.\$c
#SBATCH --output=\${ROOT_PATH}/logs/run_hinta_chr$i.\$c.out
#SBATCH --error=\${ROOT_PATH}/logs/run_hinta_chr$i.\$c.err
#SBATCH --time=23:59:00
#SBATCH --mem=20gb

python /path/to/hinta_script.py \\
  \${ROOT_PATH}/input/chunks/chr${i}.offspring.\$c.haps \\
  \${ROOT_PATH}/input/chunks/offspring.\$c.sample \\
  150 \\
  \${ROOT_PATH}/output/chunks/chr${i}.\$c.transmitted.haps \\
  \${ROOT_PATH}/output/chunks/chr${i}.\$c.transmitted.sample \\
  \${ROOT_PATH}/output/chunks/chr${i}.\$c.nontransmitted.haps \\
  \${ROOT_PATH}/output/chunks/chr${i}.\$c.nontransmitted.sample
  \" > script_hinta_chr$i.\$c.sh
  sbatch script_hinta_chr$i.\$c.sh
done" > script_hinta_chr$i.sh

  sbatch script_hinta_chr$i.sh
done
