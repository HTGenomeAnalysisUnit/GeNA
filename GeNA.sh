#!/bin/bash
set -eo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

### Parse arguments
covs=""
corr_batch="False"                                                                                                                                                               
ks=""
sampleid="id"
cleanup="Y"
connectivities_prefix=""
while [[ "$#" -gt 0 ]]
  do
    case $1 in
      -s|--single_cell_data_object) sc_object="$2"; shift;;
      -g|--genotypes) gtypes="$2"; shift;;
      -i|--sampleid) sampleid="$2"; shift;;
      -u|--connectivities) connectivities_prefix="$2"; shift;;
      -o|--output_folder) res_folder="$2"; shift;;
      -c|--covs) covs="$2"; shift;;
      -b|--corr_batch) corr_batch="$2"; shift;;
      -k|--ks) ks="$2"; shift;;
      -r|--remove_tmp) cleanup="$2"; shift;;
    esac
    shift
done

cmd_log_file="${res_folder}/commands.sh"
args_str="Genotypes: $gtypes \nSingle-cell data object: $sc_object \nSample id: $sampleid \nOutput folder: $res_folder \nCorrect for batch: $corr_batch \nCleanup: $cleanup \nCommand log: $cmd_log_file \n"
if [[ -z "$covs" ]]
then
    args_str+="Covariates: None\n"
else
    args_str+="Covariates: $covs\n"
fi
if [[ -z "$ks" ]]
then
    args_str+="Values of k: Default\n\n\n"
else
    args_str+="Values of k: $ks\n\n\n"
fi
printf "Arguments provided:\n"
printf "$args_str"

### Format inputs to PLINK
echo "============================================="
echo "[$(date)] - STEP1. Formatting input for PLINK"
echo "# STEP1. Formatting input for PLINK" > $cmd_log_file

command="python3 -u ${SCRIPT_DIR}/export_nampcs.py --res_folder ${res_folder} --sc_object_path ${sc_object} --sampleid ${sampleid}"
if [[ "$corr_batch" == "True" ]]
then
    command+=" --corr_batch 'True'"
fi
if [[ ! -z "$covs" ]]
then
    command+=" --covs $covs"
fi
if [[ ! -z "$ks" ]]
then
    command+=" --ks $ks"
fi
if [[ ! -z "$connectivities_prefix" ]]
then
    command+=" --use_custom_connectivities $connectivities_prefix"
fi

echo $command
echo $command >> $cmd_log_file
eval $command
k_max=$(awk 'NR==1{max = $1 + 0; next} {if ($1 > max) max = $1;} END {print max}' ${res_folder}/ks.csv)

### Applies PLINK to generate a test statistic reflecting the relationship between each allele and a single NAMPC
echo "============================================="
echo "[$(date)] - STEP2. Run plink association for NAMPC"
echo "\n# STEP2. Run plink association for NAMPC" >> $cmd_log_file

mkdir -p ${res_folder}/plink_per_nampc
command="plink2 --pfile ${gtypes} --pheno ${res_folder}/nampcs.csv --glm allow-no-covars --prune --out ${res_folder}/plink_per_nampc/NAM"

echo $command
echo $command >> $cmd_log_file
eval $command

command="paste"
for n_nampc in $(eval echo "{1..$k_max}")
do
    command+=" <(awk 'NR>1 {print \$12}' ${res_folder}/plink_per_nampc/NAM.PC${n_nampc}.glm.linear )"
done
command+="> ${res_folder}/t_per_nampc.txt"
echo "Gathering metrics across NAM-PCs"
echo "# Gathering metrics across NAM-PCs" >> $cmd_log_file
echo $command >> $cmd_log_file
eval $command

### Multi-nampc test
echo "============================================="
echo "[$(date)] - STEP3. multi-NAM-PC tests"
echo "\n# STEP3. multi-NAM-PC tests"  >> $cmd_log_file

command="Rscript ${SCRIPT_DIR}/joint_test.R --outfile ${res_folder}/P_k.txt \
         --chisq_per_nampc_file ${res_folder}/t_per_nampc.txt \
         --ks_file ${res_folder}/ks.csv"

echo $command
echo $command >> $cmd_log_file
eval $command

### Assemble results file
echo "============================================="
echo "[$(date)] - STEP4. Assembling GeNA results file"
echo "\n# STEP4. Assembling GeNA results file" >> $cmd_log_file

command="paste"
for i_col in $(eval echo "{1..5} 7 9 10 11") # SNP information columns
do
    command+=" <(awk '{print \$$i_col}' ${res_folder}/plink_per_nampc/NAM.PC1.glm.linear)"
done
command+=" ${res_folder}/P_k.txt" 

for n_nampc in $(eval echo "{1..$k_max}")
do
    awk_str='"BETA", "BETA_NAMPC'
    awk_str+=$(echo $n_nampc)
    awk_str+='"'
    command+=" <(awk '{print \$12}' ${res_folder}/plink_per_nampc/NAM.PC${n_nampc}.glm.linear | awk '(NR==1){gsub($awk_str, \$0);}{print;}' )"
done
command+=" > ${res_folder}/GeNA_sumstats.txt"

echo $command >> $cmd_log_file                                                                                              
eval $command

echo "============================================="
echo "[$(date)] - STEP5. Loci pruning"
echo "\n# STEP5. Loci pruning" >> $cmd_log_file

command="plink2 --pfile ${gtypes} --clump ${res_folder}/GeNA_sumstats.txt --clump-p1 1e-5 --clump-kb 500 --clump-p2 1e-3 --out ${res_folder}/GeNA_sumstats.clump"
echo $command >> $cmd_log_file                                                                                              
eval $command

echo "Extract significant loci with index SNP P < 5e-8 and at least 5 SNPs"
echo "\n# Extract significant loci with index SNP P < 5e-8 and at least 5 SNPs" >> $cmd_log_file 
command="awk 'NR == 1 || (\$5 >= 5 && \$4 <= 5e-8)' ${res_folder}/GeNA_sumstats.clump.clumps > ${res_folder}/GeNA_sumstats.clump.clumps.filtered.tsv"
echo $command >> $cmd_log_file                                                                                              
eval $command

command="cut -f3 ${res_folder}/GeNA_sumstats.clump.clumps.filtered.tsv | tail -n+2 > ${res_folder}/GeNA_sumstats.clump.clumps.filtered.indexsnps"
echo $command >> $cmd_log_file                                                                                              
eval $command

echo "============================================="
echo "[$(date)] - STEP6. Extract index SNPs genotypes"
echo "\n# STEP6. Extract index SNPs genotypes" >> $cmd_log_file

command="plink2 --pfile ${gtypes} --extract ${res_folder}/GeNA_sumstats.clump.clumps.filtered.indexsnps --export A --out ${res_folder}/GeNA_sumstats.clump.clumps.filtered.indexsnps.recodeA"
echo $command >> $cmd_log_file                                                                                              
eval $command

### Clean up
if [[ "$cleanup" == "Y" ]] 
then
  echo "Clean up temp files"
  command="rm ${res_folder}/P_k.txt"
  eval $command
  command="rm ${res_folder}/t_per_nampc.txt"
  eval $command
fi
