# intersect dot location

module load bedtools/2.28.0

loop_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/snakedots/;
cooler_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/results/comb_replicates_coolers/;
savepath=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/snakedots/comp_dots;
script_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/snakedots/comp_dots_scripts;

filename=$(ls $loop_path | grep "hg38" |  sort -r  )
#echo $filename
k=2
for i in $(seq 1 9); do
for j in $(seq $k 10);do
file1=$(echo $filename | cut -f$i -d" ");
file2=$(echo $filename | cut -f$j -d" ");
echo $file1 $file2;
#echo $file2;
f1=$(echo $file1 | cut -f1 -d".");
f2=$(echo $file2 | cut -f1 -d".");
comb_files=$f1"_vs_"$f2".txt";
echo $comb_files;
python $script_path/comp_dots.py $loop_path $cooler_path $savepath $file1 $file2
done;
k=$(echo $((k + 1)));done
