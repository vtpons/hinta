#!/bin/bash
#SBATCH --job-name=split_and_convert_to_vcf
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb
#SBATCH --nodes=1

ROOT_PATH=/path/to/project

for i in {1..22}
do
echo "#!/bin/bash
#SBATCH --job-name=process_chr$i
#SBATCH --output=${ROOT_PATH}/logs/process_files_chr$i.out
#SBATCH --error=${ROOT_PATH}/logs/process_files_chr$i.err
#SBATCH --time=20:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=10gb
#SBATCH --nodes=1

module load BCFtools

ROOT_PATH=/path/to/project

# 1. Duplicate transmitted haps/sample
echo \"\$(date): duplicating transmitted haps/sample...\" && \\

awk '{
    printf \"%s %s %s %s %s\", \$1, \$2, \$3, \$4, \$5;
    for (i=6; i<=NF; i++) {
        if (\$i == \"N\") {
            j = (i % 2 == 0) ? i-1 : i+1;
            j = (j <= NF) ? j : i;
            printf \" %s %s\", \$j, \$j;
        } else {
            printf \" %s %s\", \$i, \$i;
        }
    }
    print \"\";
}' \${ROOT_PATH}/output/chr${i}.transmitted.haps > \${ROOT_PATH}/prs_input/transmitted/chr${i}.transmitted.split.haps

awk 'NR<=2 {print \$0; next} {
    print \$0
    printf \"%s_2 %s_2\", \$1, \$2
    for (i=3; i<=NF; i++) printf \" %s\", \$i
    print \"\"
}' \${ROOT_PATH}/output/chr${i}.transmitted.sample > \${ROOT_PATH}/prs_input/transmitted/chr${i}.transmitted.split.sample && \\

# 2. Duplicate non-transmitted haps/sample
echo \"\$(date): duplicating nontransmitted haps/sample...\" && \\

awk '{
    printf \"%s %s %s %s %s\", \$1, \$2, \$3, \$4, \$5;
    for (i=6; i<NF; i+=2) {
        val1 = (\$i == \"N\") ? \$(i+1) : \$i;
        val2 = (\$(i+1) == \"N\") ? \$i : \$(i+1);
        printf \" %s %s %s %s\", val1, val1, val2, val2;
    }
    print \"\";
}' \${ROOT_PATH}/output/chr${i}.nontransmitted.haps > \${ROOT_PATH}/prs_input/nontransmitted/chr${i}.nontransmitted.split.haps

awk 'NR<=2 {print \$0; next} {
    print \$0
    printf \"%s_2 %s_2\", \$1, \$2
    for (i=3; i<=NF; i++) printf \" %s\", \$i
    print \"\"
}' \${ROOT_PATH}/output/chr${i}.nontransmitted.sample > \${ROOT_PATH}/prs_input/nontransmitted/chr${i}.nontransmitted.split.sample && \\

# 3. Convert to VCF
echo \"\$(date): converting to VCF...\" && \\

shapeit -convert --input-haps \${ROOT_PATH}/prs_input/transmitted/chr${i}.transmitted.split \\
  --output-vcf \${ROOT_PATH}/prs_input/transmitted/chr${i}.transmitted.split.vcf && \\
bgzip \${ROOT_PATH}/prs_input/transmitted/chr${i}.transmitted.split.vcf && \\
bcftools index -f \${ROOT_PATH}/prs_input/transmitted/chr${i}.transmitted.split.vcf.gz && \\

shapeit -convert --input-haps \${ROOT_PATH}/prs_input/nontransmitted/chr${i}.nontransmitted.split \\
  --output-vcf \${ROOT_PATH}/prs_input/nontransmitted/chr${i}.nontransmitted.split.vcf && \\
bgzip \${ROOT_PATH}/prs_input/nontransmitted/chr${i}.nontransmitted.split.vcf && \\
bcftools index -f \${ROOT_PATH}/prs_input/nontransmitted/chr${i}.nontransmitted.split.vcf.gz
" > script_process_chr${i}.sh

sbatch script_process_chr${i}.sh
done
