#!/bin/bash
#SBATCH --job-name=merge_hinta_output
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb
#SBATCH --nodes=1

ROOT_PATH=/path/to/project

for i in {1..22}
do
echo "#!/bin/bash
#SBATCH --job-name=merge_chr$i
#SBATCH --output=${ROOT_PATH}/logs/merge_chr$i.out
#SBATCH --error=${ROOT_PATH}/logs/merge_chr$i.err
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=5gb
#SBATCH --nodes=1

input_path=${ROOT_PATH}/output/subsets
output_path=${ROOT_PATH}/output

merge_samples() {
    local path=\$1
    local prefix=\$2
    local suffix=\$3
    local output_file=\$4

    # Take header from the first file
    head -2 \${path}/\${prefix}0.\${suffix} > \${output_file}

    for i in \$(seq 0 64); do
        awk 'FNR>2' \${path}/\${prefix}\${i}.\${suffix} >> \${output_file}
    done
}

merge_samples \"\${input_path}\" \"chr${i}.\" \"transmitted.sample\" \"\${output_path}/chr${i}.transmitted.sample\" && 
merge_samples \"\${input_path}\" \"chr${i}.\" \"nontransmitted.sample\" \"\${output_path}/chr${i}.nontransmitted.sample\" && 

merge_haps() {
    local input_pattern=\$1
    local output_file=\$2

    files=( \$(ls \${input_pattern} | sort -V) ) 

    awk '{print \$1 \" \" \$2 \" \" \$3 \" \" \$4 \" \" \$5}' \${files[0]} > \${output_file}

    for file in \"\${files[@]:0}\"; do
        paste -d \" \" \${output_file} <(awk '{for(i=6; i<NF; i++) printf \"%s \", \$i; print \$NF}' \${file}) > \${input_path}/temp_merged_$i.haps
        mv \${input_path}/temp_merged_$i.haps \${output_file}
    done
}

merge_haps \"\${input_path}/chr${i}.*.transmitted.haps\" \"\${output_path}/chr${i}.transmitted.haps\" &&
merge_haps \"\${input_path}/chr${i}.*.nontransmitted.haps\" \"\${output_path}/chr${i}.nontransmitted.haps\" 
" > script_merge_chr${i}.sh

sbatch script_merge_chr${i}.sh
done
