import cooler
import numpy as np
import pandas as pd
import glob, os
import matplotlib.pyplot as plt

#os.chdir("/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/results/coolers_library_group")
# bsub -n 2 -o subcompartments_input.log -e subcompartments_input.err -q short -W 4:00 -R rusage[mem=50000] "python subcompartments_input.py"


path="/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/results/coolers_library/"

for file in glob.glob(path+"*Rao2014-GM12878-MboI-allreps-filtered.100kb.cool"):
	c = cooler.Cooler(file)#+'::/resolutions/100000')
	name=file.split('/')[-1]
	name=name.split('.')[0]
	print(name)
	odd_chroms=c.chromnames[:-3:2]
	even_chroms=c.chromnames[1:-3:2]
	even_list=[]
	for odd_chr in odd_chroms:
		odd_list=[]
		for even_chr in even_chroms:
			chr_split=c.matrix(balance=True).fetch(even_chr,odd_chr)
			#print(chr_split.shape)
			odd_list.append(chr_split)
		odd_concat=np.concatenate(odd_list,axis=0)
		even_list.append(odd_concat)
		print("    {}".format(odd_concat.shape))
	even_odd=np.concatenate(even_list,axis=1)
	print(even_odd.shape)
	# modify that with savetxt or npz
	np.savez('/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/subcompartment_clustering/'+name+'_100kb_inter', even_odd)
	plt.clf()
	plt.imshow( np.log(even_odd),cmap='autumn')
	plt.savefig('/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_deep/subcompartment_clustering/'+name+'_100kb_inter.pdf')
