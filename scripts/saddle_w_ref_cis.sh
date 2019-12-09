 # bsub -n 1 -o saddle_w_ref_cis.log -e saddle_w_ref_cis.err -q long -W 15:00 -R rusage[mem=25000] "bash saddle_w_ref_cis.sh"

input_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/cool_matrix
output_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/cooler_snakemake/mapq30/compute_saddle_w_ref_cis_rebin
expected_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/cooler_snakemake/mapq30/expected
comp_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/cooler_snakemake/mapq30/call_compartments/U54-H1ESC4DN-FA-DSG-MNase-R1.hg38.100000.mapq_30.cis.vecs.tsv

cd $output_path
for i in $(ls $input_path | grep "mapq_30.1000.mcool" | grep "ESC-" | grep "T1") ; do \
filename=$(echo $i | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i
#compartment=$comp_path/$filename".vecs.tsv"
expected=$expected_path/$filename".hg38.100000.mapq_30.cis.expected.tsv"
echo $input
echo $comp_path
echo $expected
if [ ! -f $output_path/$filename"_cis_ESC_FA_MNase_deep.pdf" ]; then
	cooltools compute-saddle --strength --quantiles --vmin 0.5 --vmax 2 --qrange 0.02 0.98 --contact-type cis --out-prefix $output_path/$filename"_cis_ESC_FA_MNase_deep" --fig pdf $input::/resolutions/100000 $comp_path $expected;
fi;done

###

input_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/cool_matrix
output_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/cooler_snakemake/mapq30/compute_saddle_w_ref_cis_rebin
expected_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/cooler_snakemake/mapq30/expected
comp_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/cooler_snakemake/mapq30/call_compartments/U54-ESC4DN-FA-DpnII-2017524-R1-T1__hg38.hg38.100000.mapq_30.cis.vecs.tsv

cd $output_path
for i in $(ls $input_path | grep "mapq_30.1000.mcool" | grep "ESC-" | grep "T1") ; do \
filename=$(echo $i | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i
#compartment=$comp_path/$filename".vecs.tsv"
expected=$expected_path/$filename".hg38.100000.mapq_30.cis.expected.tsv"
echo $input
echo $comp_path
echo $expected
if [ ! -f $output_path/$filename"_cis_ESC_FA_DpnII_deep.pdf" ]; then
	cooltools compute-saddle --strength --quantiles --vmin 0.5 --vmax 2 --qrange 0.02 0.98  --contact-type cis --out-prefix $output_path/$filename"_cis_ESC_FA_DpnII_deep" --fig pdf $input::/resolutions/100000 $comp_path $expected;
fi;done
##########

input_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/cool_matrix
output_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/cooler_snakemake/mapq30/compute_saddle_w_ref_cis_rebin
expected_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/cooler_snakemake/mapq30/expected
comp_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/cooler_snakemake/mapq30/call_compartments/U54-ESC4DN-DSG-DpnII-20190530-R1-T1__hg38.hg38.100000.mapq_30.cis.vecs.tsv

cd $output_path
for i in $(ls $input_path | grep "mapq_30.1000.mcool" | grep "ESC-" | grep "T1") ; do \
filename=$(echo $i | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i
#compartment=$comp_path/$filename".vecs.tsv"
expected=$expected_path/$filename".hg38.100000.mapq_30.cis.expected.tsv"
echo $input
echo $comp_path
echo $expected
if [ ! -f $output_path/$filename"_cis_ESC_DSG_DpnII_deep.pdf" ]; then
	cooltools compute-saddle --strength --quantiles --vmin 0.5 --vmax 2 --qrange 0.02 0.98  --contact-type cis --out-prefix $output_path/$filename"_cis_ESC_DSG_DpnII_deep" --fig pdf $input::/resolutions/100000 $comp_path $expected;
fi;done



##

input_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/cool_matrix
output_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/cooler_snakemake/mapq30/compute_saddle_w_ref_cis_rebin
expected_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/cooler_snakemake/mapq30/expected
comp_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/cooler_snakemake/mapq30/call_compartments/U54-END4DN-FA-DSG-MNase-R1.hg38.100000.mapq_30.cis.vecs.tsv

cd $output_path
for i in $(ls $input_path | grep "mapq_30.1000.mcool" | grep "END-" | grep "T1") ; do \
filename=$(echo $i | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i
#compartment=$comp_path/$filename".vecs.tsv"
expected=$expected_path/$filename".hg38.100000.mapq_30.cis.expected.tsv"
echo $input
echo $comp_path
echo $expected
if [ ! -f $output_path/$filename"_cis_END_FA_MNase_deep.pdf" ]; then
	cooltools compute-saddle --strength --quantiles --vmin 0.5 --vmax 2 --qrange 0.02 0.98  --contact-type cis --out-prefix $output_path/$filename"_cis_END_FA_MNase_deep" --fig pdf $input::/resolutions/100000 $comp_path $expected;
fi;done

####
input_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/cool_matrix
output_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/cooler_snakemake/mapq30/compute_saddle_w_ref_cis_rebin
expected_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/cooler_snakemake/mapq30/expected
comp_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/cooler_snakemake/mapq30/call_compartments/U54-END-FA-DpnII-030917__hg38.hg38.100000.mapq_30.cis.vecs.tsv

cd $output_path
for i in $(ls $input_path | grep "mapq_30.1000.mcool" | grep "END-" | grep "T1") ; do \
filename=$(echo $i | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i
#compartment=$comp_path/$filename".vecs.tsv"
expected=$expected_path/$filename".hg38.100000.mapq_30.cis.expected.tsv"
echo $input
echo $comp_path
echo $expected
if [ ! -f $output_path/$filename"_cis_END_FA_DpnII_deep.pdf" ]; then
	cooltools compute-saddle --strength --quantiles --vmin 0.5 --vmax 2 --qrange 0.02 0.98  --contact-type cis --out-prefix $output_path/$filename"_cis_END_FA_DpnII_deep" --fig pdf $input::/resolutions/100000 $comp_path $expected;
fi;done

####


input_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/cool_matrix
output_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/cooler_snakemake/mapq30/compute_saddle_w_ref_cis_rebin
expected_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/cooler_snakemake/mapq30/expected
comp_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/cooler_snakemake/mapq30/call_compartments/U54-HFFc64DN-FA-DSG-MNase-R1.hg38.100000.mapq_30.cis.vecs.tsv


cd $output_path
for i in $(ls $input_path | grep "mapq_30.1000.mcool" | grep "HFF-") ; do \
filename=$(echo $i | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i
#compartment=$comp_path/$filename".vecs.tsv"
expected=$expected_path/$filename".hg38.100000.mapq_30.cis.expected.tsv"
echo $input
echo $comp_path
echo $expected
if [ ! -f $output_path/$filename"_cis_HFFc6_FA_MNase_deep.pdf" ]; then
	cooltools compute-saddle --strength --quantiles --vmin 0.5 --vmax 2 --qrange 0.02 0.98  --contact-type cis --out-prefix $output_path/$filename"_cis_HFFc6_FA_MNase_deep" --fig pdf $input::/resolutions/100000 $comp_path $expected;
fi;done


###


input_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/cool_matrix
output_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/cooler_snakemake/mapq30/compute_saddle_w_ref_cis_rebin
expected_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/cooler_snakemake/mapq30/expected
comp_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/cooler_snakemake/mapq30/call_compartments/U54-HFFc6-p17-FA-DpnII-20170327__hg38.hg38.100000.mapq_30.cis.vecs.tsv


cd $output_path
for i in $(ls $input_path | grep "mapq_30.1000.mcool" | grep "HFF-") ; do \
filename=$(echo $i | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i
#compartment=$comp_path/$filename".vecs.tsv"
expected=$expected_path/$filename".hg38.100000.mapq_30.cis.expected.tsv"
echo $input
echo $comp_path
echo $expected
if [ ! -f $output_path/$filename"_cis_HFFc6_FA_DpnII_deep.pdf" ]; then
	cooltools compute-saddle --strength --quantiles --vmin 0.5 --vmax 2 --qrange 0.02 0.98 --contact-type cis --out-prefix $output_path/$filename"_cis_HFFc6_FA_DpnII_deep" --fig pdf $input::/resolutions/100000 $comp_path $expected;
fi;done


# ###

input_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/cool_matrix
output_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/cooler_snakemake/mapq30/compute_saddle_w_ref_cis_rebin
expected_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/cooler_snakemake/mapq30/expected
comp_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/cooler_snakemake/mapq30/call_compartments/U54-HFFc6-DSG-DpnII-20180319-R1-T1__hg38.hg38.100000.mapq_30.cis.vecs.tsv

cd $output_path
for i in $(ls $input_path | grep "mapq_30.1000.mcool" | grep "HFF-") ; do \
filename=$(echo $i | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i
#compartment=$comp_path/$filename".vecs.tsv"
expected=$expected_path/$filename".hg38.100000.mapq_30.cis.expected.tsv"
echo $input
echo $comp_path
echo $expected
if [ ! -f $output_path/$filename"_cis_HFFc6_DSG_DpnII_deep.pdf" ]; then
	cooltools compute-saddle --strength --quantiles --vmin 0.5 --vmax 2 --qrange 0.02 0.98 --contact-type cis --out-prefix $output_path/$filename"_cis_HFFc6_DSG_DpnII_deep" --fig pdf $input::/resolutions/100000 $comp_path $expected;
fi;done

##

input_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/cool_matrix
output_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/cooler_snakemake/mapq30/compute_saddle_w_ref_cis_rebin
expected_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/cooler_snakemake/mapq30/expected
comp_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/cooler_snakemake/mapq30/call_compartments/U54-HFFc6-DSG-DdeI-20180319-R1-T1__hg38.hg38.100000.mapq_30.cis.vecs.tsv

cd $output_path
for i in $(ls $input_path | grep "mapq_30.1000.mcool" | grep "HFF-") ; do \
filename=$(echo $i | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i
#compartment=$comp_path/$filename".vecs.tsv"
expected=$expected_path/$filename".hg38.100000.mapq_30.cis.expected.tsv"
echo $input
echo $comp_path
echo $expected
if [ ! -f $output_path/$filename"_cis_HFFc6_DSG_DdeI_deep.pdf" ]; then
	cooltools compute-saddle --strength --quantiles --vmin 0.5 --vmax 2 --qrange 0.02 0.98 --contact-type cis --out-prefix $output_path/$filename"_cis_HFFc6_DSG_DdeI_deep" --fig pdf $input::/resolutions/100000 $comp_path $expected;
fi;done

##

input_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/cool_matrix
output_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/cooler_snakemake/mapq30/compute_saddle_w_ref_cis_rebin
expected_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/cooler_snakemake/mapq30/expected
comp_path=/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/cooler_snakemake/mapq30/call_compartments/U54-HFFc6-DSG-DdeI-DpnII-20190711-R1-T1__hg38.hg38.100000.mapq_30.cis.vecs.tsv

cd $output_path
for i in $(ls $input_path | grep "mapq_30.1000.mcool" | grep "HFF-") ; do \
filename=$(echo $i | cut -f1 -d".");
#echo $i;done
echo $filename
input=$input_path/$i
#compartment=$comp_path/$filename".vecs.tsv"
expected=$expected_path/$filename".hg38.100000.mapq_30.cis.expected.tsv"
echo $input
echo $comp_path
echo $expected
if [ ! -f $output_path/$filename"_cis_HFFc6_DSG_DdeI_DpnII_deep.pdf" ]; then
	cooltools compute-saddle --strength --quantiles --vmin 0.5 --vmax 2 --qrange 0.02 0.98 --contact-type cis --out-prefix $output_path/$filename"_cis_HFFc6_DSG_DdeI_DpnII_deep" --fig pdf $input::/resolutions/100000 $comp_path $expected;
fi;done
