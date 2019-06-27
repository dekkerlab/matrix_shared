#bsub -n 1 -o submit_compute_comp_per_arm.log -e submit_compute_comp_per_arm.err -q long -W 20:00 -R rusage[mem=100000] "bash submit_compute_comp_per_arm.sh"
module load bedtools/2.28.0

input_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/cool_matrix
output_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/cool_matrix
cd $output_path
for i in $(ls $input_path | grep "mapq_30.1000.mcool") ; do \
filename=$(echo $i | cut -f1 -d"_");
echo $filename;
if [ ! -f $filename".cis.250000.bw" ]; then
python /nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/cooler_snakemake/comp_per_arm/compute_comp_per_arm.py \
$input_path/$i 250000 ;fi;done
