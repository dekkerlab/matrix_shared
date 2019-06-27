
## 
list_stats = list.files(path = "/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/results/pairs_library", pattern ="dedup.stats")
stat_file = c()
for (i in 1:length(list_stats))
{
	
	f = read.table(paste("/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/results/pairs_library/",list_stats[i],sep=""),header=F, sep="\t")
	rownames(f) <- f[,1]

	total_reads=f["total",2]
	total_unmapped=f["total_unmapped",2]
	total_mapped=f["total_mapped",2]
	total_nodups=f["total_nodups",2]
	cis=f["cis",2]
	cis1kb=f["cis_1kb+",2]
	trans=f["trans",2]
	mapped_ratio=100*(total_nodups/total_reads)
	cis1kb_ratio=100*(cis1kb/total_nodups)
	cis_ratio=100*(cis/total_nodups)
	trans_ratio=100*(trans/total_nodups)
	one_side_mapped=f[3,2]
	
	row = cbind(list_stats[i],total_reads,total_nodups,total_unmapped,mapped_ratio,cis,cis_ratio,cis1kb_ratio,trans,trans_ratio,one_side_mapped)
	stat_file=rbind(stat_file,row)

}
write.table(stat_file,"Stats_all_samples.txt",sep="\t", row.names=F)


### Distribution of cis after diagonal removal 

pattern='ESC'

library(RColorBrewer)

list_stats = intersect(list.files(path = "/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/results/pairs_library", pattern =pattern),list.files(path = "/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/results/pairs_library", pattern ="dedup.stats"))
stat_file = c()
l = length(list_stats)
print(list_stats)

pdf(paste("Cis_percent_with_distance",pattern,".pdf",sep=""))
color=rainbow(l)
k=3
l=3
m=3
FA_col=brewer.pal(n = 6, name = 'Blues')
DSG_col=brewer.pal(n = 6, name = 'Greens')
EGS_col=brewer.pal(n = 6, name = 'Reds')

for (i in 1:length(list_stats))
{
	f = read.table(paste("/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/results/pairs_library/",list_stats[i],sep=""),header=F, sep="\t")
	rownames(f) <- f[,1]
	total_nodups=f["total_nodups",2]
	cis=f["cis",2]
	cis1kb=f["cis_1kb+",2]
	cis2kb=f["cis_2kb+",2]
	cis4kb=f["cis_4kb+",2]
	cis10kb=f["cis_10kb+",2]
	cis20kb=f["cis_20kb+",2]
	cis40kb=f["cis_40kb+",2]

	cis_ratio=100*(cis/total_nodups)
	cis1kb_ratio=100*(cis1kb/total_nodups)
	cis2kb_ratio=100*(cis2kb/total_nodups)
	cis4kb_ratio=100*(cis4kb/total_nodups)
	cis10kb_ratio=100*(cis10kb/total_nodups)
	cis20kb_ratio=100*(cis20kb/total_nodups)
	cis40kb_ratio=100*(cis40kb/total_nodups)
	
	xkb <- c(0,1,2,3,4,5,6)
	xkb_names <- c("0kb","1kb","2kb","4kb","10kb","20kb","40kb")
	row = cbind(list_stats[i],total_nodups,cis,cis_ratio,cis1kb_ratio,cis2kb_ratio,cis4kb_ratio,cis10kb_ratio,cis20kb_ratio,cis40kb_ratio)
	stat_file=rbind(stat_file,row)
	#print(dim(stat_file[1]))
	if(dim(stat_file)[1]==1)
	{
		if (grepl("FA", list_stats[i]))
		{
		plot(xkb,as.numeric(stat_file[i,4:10]),ylim=range(0,100),cex=0.5,pch=16,col=FA_col[k],ylab="Cis Percent",xaxt="n",xlab="kb")
		lines(xkb,as.numeric(stat_file[i,4:10]),col=FA_col[k])
		text(x=xkb,y=-5,labels=xkb_names, srt=45, adj=1, xpd=TRUE)
		k=k+1
		print("FA")
		}

        if (grepl("DSG", list_stats[i]))
		{
		plot(xkb,as.numeric(stat_file[i,4:10]),ylim=range(0,100),cex=0.5,pch=16,col=DSG_col[l],ylab="Cis Percent",xaxt="n",xlab="kb")
		lines(xkb,as.numeric(stat_file[i,4:10]),col=DSG_col[l])
		text(x=xkb,y=-5,labels=xkb_names, srt=45, adj=1, xpd=TRUE)
		l=l+1
		print("DSG")
		}

		if (grepl("EGS", list_stats[i]))
		{
		plot(xkb,as.numeric(stat_file[i,4:10]),ylim=range(0,100),cex=0.5,pch=16,col=EGS_col[m],ylab="Cis Percent",xaxt="n",xlab="kb")
		lines(xkb,as.numeric(stat_file[i,4:10]),col=EGS_col[m])
		text(x=xkb,y=-5,labels=xkb_names, srt=45, adj=1, xpd=TRUE)
		m=m+1
		print("EGS")
		}
	}
	if(dim(stat_file)[1]>1)
	{
		if (grepl("FA", list_stats[i]))
		{
		points(xkb,as.numeric(stat_file[i,4:10]),pch=16,cex=0.5,col=FA_col[k])
		lines(xkb,as.numeric(stat_file[i,4:10]),col=FA_col[k])
		k=k+1
		}

        if (grepl("DSG", list_stats[i]))
		{
		points(xkb,as.numeric(stat_file[i,4:10]),pch=16,cex=0.5,col=DSG_col[l])
		lines(xkb,as.numeric(stat_file[i,4:10]),col=DSG_col[l])
		l=l+1
		}

		if (grepl("EGS", list_stats[i]))
		{
		points(xkb,as.numeric(stat_file[i,4:10]),pch=16,cex=0.5,col=EGS_col)
		lines(xkb,as.numeric(stat_file[i,4:10]),col=EGS_col[m])
		m=m+1
		}

	}

}
names=stat_file[,1]
names_new=gsub("-2.*","",names)
names_new=gsub("-0.*","",names_new)
legend("topright", names_new,fill=c(DSG_col[3:6],EGS_col[3:6],FA_col[3:6]),cex=0.5,border=c(DSG_col[3:6],EGS_col[3:6],FA_col[3:6]))
dev.off()



### Distribution of cis after diagonal removal ##  recalculate the cis percentage


pattern="ESC"
library(RColorBrewer)

list_stats = intersect(list.files(path = "/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/results/pairs_library", pattern =pattern),list.files(path = "/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/results/pairs_library", pattern ="dedup.stats"))
stat_file = c()
l = length(list_stats)
print(list_stats)

pdf(paste("Cis_percent_with_distance",pattern,"_newcis.pdf",sep=""))
color=rainbow(l)
k=3
l=3
m=3
FA_col=brewer.pal(n = 6, name = 'Blues')
DSG_col=brewer.pal(n = 6, name = 'Greens')
EGS_col=brewer.pal(n = 6, name = 'Reds')

for (i in 1:length(list_stats))
{
	f = read.table(paste("/nl/umw_job_dekker/users/ba69w/HiC_Analysis/U54_matrix/results/pairs_library/",list_stats[i],sep=""),header=F, sep="\t")
	rownames(f) <- f[,1]
	total_nodups=f["total_nodups",2]
	cis=f["cis",2]
	trans=f["trans",2]
	cis1kb=f["cis_1kb+",2]
	cis2kb=f["cis_2kb+",2]
	cis4kb=f["cis_4kb+",2]
	cis10kb=f["cis_10kb+",2]
	cis20kb=f["cis_20kb+",2]
	cis40kb=f["cis_40kb+",2]

	cis_ratio=100*(cis/total_nodups)
	cis1kb_ratio=100*(cis1kb/(cis1kb+trans))
	cis2kb_ratio=100*(cis2kb/(cis2kb+trans))
	cis4kb_ratio=100*(cis4kb/(cis4kb+trans))
	cis10kb_ratio=100*(cis10kb/(cis10kb+trans))
	cis20kb_ratio=100*(cis20kb/(cis20kb+trans))
	cis40kb_ratio=100*(cis40kb/(cis40kb+trans))
	print(cis+trans)
	print(total_nodups)
	xkb <- c(0,1,2,3,4,5,6)
	xkb_names <- c("0kb","1kb","2kb","4kb","10kb","20kb","40kb")
	row = cbind(list_stats[i],total_nodups,cis,cis_ratio,cis1kb_ratio,cis2kb_ratio,cis4kb_ratio,cis10kb_ratio,cis20kb_ratio,cis40kb_ratio)
	stat_file=rbind(stat_file,row)
	#print(dim(stat_file[1]))
	if(dim(stat_file)[1]==1)
	{
		if (grepl("FA", list_stats[i]))
		{
		plot(xkb,as.numeric(stat_file[i,4:10]),ylim=range(0,100),cex=0.5,pch=16,col=FA_col[k],ylab="Cis Percent",xaxt="n",xlab="kb")
		lines(xkb,as.numeric(stat_file[i,4:10]),col=FA_col[k])
		text(x=xkb,y=-5,labels=xkb_names, srt=45, adj=1, xpd=TRUE)
		k=k+1
		print("FA")
		}

        if (grepl("DSG", list_stats[i]))
		{
		plot(xkb,as.numeric(stat_file[i,4:10]),ylim=range(0,100),cex=0.5,pch=16,col=DSG_col[l],ylab="Cis Percent",xaxt="n",xlab="kb")
		lines(xkb,as.numeric(stat_file[i,4:10]),col=DSG_col[l])
		text(x=xkb,y=-5,labels=xkb_names, srt=45, adj=1, xpd=TRUE)
		l=l+1
		print("DSG")
		}

		if (grepl("EGS", list_stats[i]))
		{
		plot(xkb,as.numeric(stat_file[i,4:10]),ylim=range(0,100),cex=0.5,pch=16,col=EGS_col[m],ylab="Cis Percent",xaxt="n",xlab="kb")
		lines(xkb,as.numeric(stat_file[i,4:10]),col=EGS_col[m])
		text(x=xkb,y=-5,labels=xkb_names, srt=45, adj=1, xpd=TRUE)
		m=m+1
		print("EGS")
		}
	}
	if(dim(stat_file)[1]>1)
	{
		if (grepl("FA", list_stats[i]))
		{
		points(xkb,as.numeric(stat_file[i,4:10]),pch=16,cex=0.5,col=FA_col[k])
		lines(xkb,as.numeric(stat_file[i,4:10]),col=FA_col[k])
		k=k+1
		}

        if (grepl("DSG", list_stats[i]))
		{
		points(xkb,as.numeric(stat_file[i,4:10]),pch=16,cex=0.5,col=DSG_col[l])
		lines(xkb,as.numeric(stat_file[i,4:10]),col=DSG_col[l])
		l=l+1
		}

		if (grepl("EGS", list_stats[i]))
		{
		points(xkb,as.numeric(stat_file[i,4:10]),pch=16,cex=0.5,col=EGS_col)
		lines(xkb,as.numeric(stat_file[i,4:10]),col=EGS_col[m])
		m=m+1
		}

	}

}
names=stat_file[,1]
names_new=gsub("-2.*","",names)
names_new=gsub("-0.*","",names_new)
legend("topright", names_new,fill=c(DSG_col[3:6],EGS_col[3:6],FA_col[3:6]),cex=0.5,border=c(DSG_col[3:6],EGS_col[3:6],FA_col[3:6]))
dev.off()



