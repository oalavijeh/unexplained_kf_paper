## PCA Case/Control Matching
## A. P. Levine // a.levine[at]ucl.ac.uk
## July 2020

#setwd("~/re_gecip/renal/oalavijeh/projects/stones/")

#Find matched controls
ac = function(x){as.character(x)}
an = function(x){as.numeric(as.character(x))}

#Load weightings
ev = an(read.table("~/re_gecip/renal/oalavijeh/projects/phenotypes/inputs_phenotypes/plink.eigenval",header=F)[,1])

#Load PC data
d=read.table("~/re_gecip/renal/oalavijeh/projects/phenotypes/inputs_phenotypes/plink.eigenvec", header=F)
names(d)=c("FID","IID",paste("PC",1:10,sep=""))

# Load ped file
a <- read.table("~/re_gecip/renal/oalavijeh/projects/phenotypes/inputs_phenotypes/final_unrelated_cystic_rd_controls.fam", header=F)
names(a) <- c("FID", "IID", "PATID", "MATID", "SEX", "COHORT")
b <- merge(a,d,by="FID", all=TRUE)
b <- b[,-c(3,4,5,7)]
d<-b[,c(1,2,4,5,6,7,8,9,10,11,12,13,3)]
names(d)=c("FID","IID",paste("PC",1:10,sep=""),"COHORT")

cases=d[which(d$COHORT==2),]
controls=d[which(d$COHORT==1),]

pdf("PCA.pdf",height=10,width=14)
par(mfrow=c(2,2))
plot(controls$PC1,controls$PC2,xlab="Principal Component 1",ylab="Principal Component 2",col=8,cex.lab=1.5)
points(cases$PC1,cases$PC2,col=2,pch=16,cex=1.2)
legend(x="bottomright",legend=c("Controls","Cases"),col=c(8,2),pch=c(1,16),cex=1.2)
plot(controls$PC3,controls$PC4,xlab="Principal Component 3",ylab="Principal Component 4",col=8,cex.lab=1.5)
points(cases$PC3,cases$PC4,col=2,pch=16,cex=1.2)
legend(x="bottomright",legend=c("Controls","Cases"),col=c(8,2),pch=c(1,16),cex=1.2)
plot(controls$PC5,controls$PC6,xlab="Principal Component 5",ylab="Principal Component 6",col=8,cex.lab=1.5)
points(cases$PC5,cases$PC6,col=2,pch=16,cex=1.2)
legend(x="topright",legend=c("Controls","Cases"),col=c(8,2),pch=c(1,16),cex=1.2)
plot(controls$PC7,controls$PC8,xlab="Principal Component 7",ylab="Principal Component 8",col=8,cex.lab=1.5)
points(cases$PC7,cases$PC8,col=2,pch=16,cex=1.2)
legend(x="topright",legend=c("Controls","Cases"),col=c(8,2),pch=c(1,16),cex=1.2)
dev.off()

cases_use=cases
controls_use=controls
store=list()
prog=0
for (i in 1:length(cases_use[,1])){
	this_prog=round(i/length(cases_use[,1])*100)
	if (this_prog>prog){prog=this_prog;print(this_prog)}
	out=c()
	for (j in 1:length(controls_use[,1])){
		this=c()
		k.count=0
		for (k in paste("PC",1:10,sep="")){
			k.count=k.count+1
			#With weighting:
			this=c(this,ev[k.count]*(an(cases_use[i,k])-an(controls_use[j,k]))^2)

			#No weighting:
			#this=c(this,(an(cases_use[i,k])-an(controls_use[j,k]))^2)
		}
		out=c(out,sum(this))
	}

	x=as.data.frame(cbind(ac(controls_use$IID),out))
	x.o=x[order(an(x$out)),]
	store[[i]] = x.o
}

save(store,file="Matched_controls_distances_PC10.RData")

#If loading from saved
load(file="Matched_controls_distances_PC10.RData")

temp=c()
count=0
for (i in paste("PC",1:10,sep="")){
	#With weighting:
	count=count+1
	temp=c(temp,ev[count]*(min(c(cases_use[,i],controls[,i]))-max(c(cases_use[,i],controls[,i])))^2)

	#No weighting:
	#temp=c(temp,(min(c(cases_use[,i],controls[,i]))-max(c(cases_use[,i],controls[,i])))^2)
}


v1=0.002 #distance from range to consider
t=v1*sum(temp)

x=2 #number of controls that must be found per case within distance t for case to be included


found=c()
remove_cases=c()
for (i in 1:length(cases_use[,1])){
	x.o = store[[i]]

	#Select the closest x controls:
	#this_found = ac(x.o[,1])[1:x]
	#or...

	#Select the x controls within a distance <=t
	this_found=ac(x.o[which(an(x.o[,2])<=t),1])

	if (length(this_found)<=x){remove_cases=c(remove_cases,i)}
	if (length(this_found)>x){found=c(found,this_found)}
	found=c(found,this_found)
}
found=unique(found)

#If using 'this_found = ac(x.o[,1])[1:x]' then:
#remove_cases=""

pdf("plot_ancestry_matched_t0.002_x2.pdf",height=10,width=14)
par(mfrow=c(2,2))
plot(controls$PC1,controls$PC2,xlab="Principal Component 1",ylab="Principal Component 2",col=8,cex.lab=1.5)
points(controls$PC1[which(controls$IID %in% found)],controls$PC2[which(controls$IID %in% found)],col=3,pch=16)
points(cases_use$PC1,cases_use$PC2,col=2,pch=16)
points(cases_use$PC1[remove_cases],cases_use$PC2[remove_cases],col=5,pch=16)
#text(cases_use$PC1[remove_cases],cases_use$PC2[remove_cases],labels=remove_cases,pch=16)
plot(controls$PC3,controls$PC4,xlab="Principal Component 3",ylab="Principal Component 4",col=8,cex.lab=1.5)
points(controls$PC3[which(controls$IID %in% found)],controls$PC4[which(controls$IID %in% found)],col=3,pch=16)
points(cases_use$PC3,cases_use$PC4,col=2,pch=16)
points(cases_use$PC3[remove_cases],cases_use$PC4[remove_cases],col=5,pch=16)
#text(cases_use$PC3[remove_cases],cases_use$PC4[remove_cases],labels=remove_cases,pch=16)
plot(controls$PC5,controls$PC6,xlab="Principal Component 5",ylab="Principal Component 6",col=8,cex.lab=1.5)
points(controls$PC5[which(controls$IID %in% found)],controls$PC6[which(controls$IID %in% found)],col=3,pch=16)
points(cases_use$PC5,cases_use$PC6,col=2,pch=16)
points(cases_use$PC5[remove_cases],cases_use$PC6[remove_cases],col=5,pch=16)
#text(cases_use$PC5[remove_cases],cases_use$PC6[remove_cases],labels=remove_cases,pch=16)
plot(controls$PC7,controls$PC8,xlab="Principal Component 7",ylab="Principal Component 8",col=8,cex.lab=1.5)
points(controls$PC7[which(controls$IID %in% found)],controls$PC8[which(controls$IID %in% found)],col=3,pch=16)
points(cases_use$PC7,cases_use$PC8,col=2,pch=16)
points(cases_use$PC7[remove_cases],cases_use$PC8[remove_cases],col=5,pch=16)
#text(cases_use$PC7[remove_cases],cases_use$PC8[remove_cases],labels=remove_cases,pch=16)
print(length(found))
print(length(remove_cases))
dev.off()

#Output table
write.table(cbind(c(ac(cases_use$IID),found),c(ac(cases_use$IID),found)),file="Matched_v15_rd_controls.txt",quote=F,row.names=F,col.names=F,sep="\t")
write.table(ac(cases_use$IID), "ancestry_matched_cystic_rd_cases.txt", sep="\t", row.names=F, col.names =F, quote=F)
write.table(cases_use$IID, "all_cystic_cases_rd.txt", sep="\t", row.names=F, col.names =F, quote=F)
