##########Rarefaction#############
##########Rarefaction#############
##########Rarefaction#############
##########Rarefaction#############

dataRaf<-read.table("Rarefaction/curacao_16S_v1v2_zotus_otutab_num_tax_20190307.txt", sep="\t",header=TRUE)
summary(dataRaf)
dim(dataRaf)

###calculate sample number

###CR1#####
dataRaf<-read.table("Rarefaction/curacao_16S_v1v2_zotus_otutab_num_tax_20190307.txt", sep="\t",header=TRUE)
sampleNumCR1<-sum(subset(dataRaf, dataRaf[,2]>=2)[,2])
dataRaf<-subset(dataRaf, dataRaf[,2]>=2)
CR1divi<-data.frame("OTU"=dataRaf[,1], "amiplicon"=dataRaf[,2], "Pi"=(dataRaf[,2]/sampleNumCR1))
summary(CR1divi)
nlevels(droplevels(CR1divi[,1]))
CR1diversity<-data.frame(CR1divi, "Piln"=(CR1divi$Pi*log(CR1divi$Pi)))
summary(CR1diversity)
nnn<-colSums(CR1diversity[,2:4])
Hmax<-log(nlevels(droplevels(CR1divi[,1])))
shannon<-abs(colSums(CR1diversity[,2:4])[3])
Evenness<-shannon/Hmax

dataRaf<-read.table("Rarefaction/curacao_16S_v1v2_zotus_otutab_num_tax_20190307.txt", sep="\t",header=TRUE)
sampleNumCR2<-sum(subset(dataRaf, dataRaf[,3]>=2)[,3])
dataRaf<-subset(dataRaf, dataRaf[,3]>=2)
CR1divi<-data.frame("OTU"=dataRaf[,1], "amiplicon"=dataRaf[,3], "Pi"=(dataRaf[,3]/sampleNumCR2))
summary(CR1divi)
nlevels(droplevels(CR1divi[,1]))
CR1diversity<-data.frame(CR1divi, "Piln"=(CR1divi$Pi*log(CR1divi$Pi)))
summary(CR1diversity)
nnn<-colSums(CR1diversity[,2:4])
Hmax<-log(nlevels(droplevels(CR1divi[,1])))
shannon<-abs(colSums(CR1diversity[,2:4])[3])
Evenness<-shannon/Hmax

dataRaf<-read.table("Rarefaction/curacao_16S_v1v2_zotus_otutab_num_tax_20190307.txt", sep="\t",header=TRUE)
sampleNumCR3<-sum(subset(dataRaf, dataRaf[,5]>=2)[,5])
dataRaf<-subset(dataRaf, dataRaf[,5]>=2)
CR1divi<-data.frame("OTU"=dataRaf[,1], "amiplicon"=dataRaf[,5], "Pi"=(dataRaf[,5]/sampleNumCR3))
summary(CR1divi)
nlevels(droplevels(CR1divi[,1]))
CR1diversity<-data.frame(CR1divi, "Piln"=(CR1divi$Pi*log(CR1divi$Pi)))
summary(CR1diversity)
nnn<-colSums(CR1diversity[,2:4])
Hmax<-log(nlevels(droplevels(CR1divi[,1])))
shannon<-abs(colSums(CR1diversity[,2:4])[3])
Evenness<-shannon/Hmax

dataRaf<-read.table("Rarefaction/curacao_16S_v1v2_zotus_otutab_num_tax_20190307.txt", sep="\t",header=TRUE)
sampleNumCR4<-sum(subset(dataRaf, dataRaf[,6]>=2)[,6])
dataRaf<-subset(dataRaf, dataRaf[,6]>=2)
CR1divi<-data.frame("OTU"=dataRaf[,1], "amiplicon"=dataRaf[,6], "Pi"=(dataRaf[,6]/sampleNumCR4))
summary(CR1divi)
nlevels(droplevels(CR1divi[,1]))
CR1diversity<-data.frame(CR1divi, "Piln"=(CR1divi$Pi*log(CR1divi$Pi)))
summary(CR1diversity)
nnn<-colSums(CR1diversity[,2:4])
Hmax<-log(nlevels(droplevels(CR1divi[,1])))
shannon<-abs(colSums(CR1diversity[,2:4])[3])
Evenness<-shannon/Hmax




###CR1#####
###CR1#####
dataRaf<-read.table("Rarefaction/curacao_16S_v1v2_zotus_otutab_num_tax_20190307.txt", sep="\t",header=TRUE)

###Create lists
CR1list<-list()
for (I in 1:dim(dataRaf)[[1]])
{
	OTU<-as.character(dataRaf[I,1])
	Numb<-dataRaf[I,2]
	if(Numb>1)
	{
		for(G in 1:Numb)
		{
			CR1list<-c(CR1list, OTU)	
		}
	}
}

####Rarefaction####
write.table(t(c("sample", "OTUs")), file="Rarefaction/CR1results.csv", sep=",", row.names=F, col.names=F)
E<-1
while(E < sampleNumCR1)
{    listed<-list()
	for (G in 1:E)
	{
		listed[G]<-CR1list[sample(1:sampleNumCR1, 1, replace=T)]
	}
	cc<-data.frame(do.call(rbind, listed))
	Value<-nlevels(cc[,1])
	goodies<-c(E, Value)
	write.table(t(goodies), "Rarefaction/CR1results.csv", sep=",", append=TRUE, row.names=F, col.names=F)
	E<-E+round(E^(1/2), digits=0)
}

###CR2#####
###CR2#####
###CR2#####

###Create lists
CR1list<-list()
for (I in 1:dim(dataRaf)[[1]])
{
	OTU<-as.character(dataRaf[I,1])
	Numb<-dataRaf[I,3]
	if(Numb>1)
	{
		for(G in 1:Numb)
		{
			CR1list<-c(CR1list, OTU)	
		}
	}
}

####Rarefaction####
write.table(t(c("sample", "OTUs")), file="Rarefaction/CR2results.csv", sep=",", row.names=F, col.names=F)
E<-1
while(E < sampleNumCR2)
{    listed<-list()
	for (G in 1:E)
	{
		listed[G]<-CR1list[sample(1:sampleNumCR2, 1, replace=T)]
	}
	cc<-data.frame(do.call(rbind, listed))
	Value<-nlevels(cc[,1])
	goodies<-c(E, Value)
	write.table(t(goodies), "Rarefaction/CR2results.csv", sep=",", append=TRUE, row.names=F, col.names=F)
	E<-E+round(E^(1/2), digits=0)
}

###CR3#####
###CR3#####
###CR3#####

###Create lists
CR1list<-list()
for (I in 1:dim(dataRaf)[[1]])
{
	OTU<-as.character(dataRaf[I,1])
	Numb<-dataRaf[I,5]
	if(Numb>1)
	{
		for(G in 1:Numb)
		{
			CR1list<-c(CR1list, OTU)	
		}
	}
}

####Rarefaction####
write.table(t(c("sample", "OTUs")), file="Rarefaction/CR3results.csv", sep=",", row.names=F, col.names=F)
E<-1
while(E < sampleNumCR3)
{    listed<-list()
	for (G in 1:E)
	{
		listed[G]<-CR1list[sample(1:sampleNumCR3, 1, replace=T)]
	}
	cc<-data.frame(do.call(rbind, listed))
	Value<-nlevels(cc[,1])
	goodies<-c(E, Value)
	write.table(t(goodies), "Rarefaction/CR3results.csv", sep=",", append=TRUE, row.names=F, col.names=F)
	E<-E+round(E^(1/2), digits=0)
}

###CR1#####
###CR1#####
###CR1#####

###Create lists
CR1list<-list()
for (I in 1:dim(dataRaf)[[1]])
{
	OTU<-as.character(dataRaf[I,1])
	Numb<-dataRaf[I,6]
	if(Numb>1)
	{
		for(G in 1:Numb)
		{
			CR1list<-c(CR1list, OTU)	
		}
	}
}

dim(CR1list)
####Rarefaction####
write.table(t(c("sample", "OTUs")), file="Rarefaction/CR4results.csv", sep=",", row.names=F, col.names=F)
E<-1
while(E < sampleNumCR4)
{   listed<-list()
	for (G in 1:E)
	{
		listed[G]<-CR1list[sample(1:sampleNumCR4, 1, replace=T)]
	}
	cc<-data.frame(do.call(rbind, listed))
	Value<-nlevels(cc[,1])
	goodies<-c(E, Value)
	write.table(t(goodies), "Rarefaction/CR4results.csv", sep=",", append=TRUE, row.names=F, col.names=F)
	E<-E+round(E^(1/2), digits=0)
}


RARE1<-read.table("Rarefaction/CR1results.csv", sep=",",header=TRUE)
RARE2<-read.table("Rarefaction/CR2results.csv", sep=",",header=TRUE)
RARE3<-read.table("Rarefaction/CR3results.csv", sep=",",header=TRUE)
RARE4<-read.table("Rarefaction/CR4results.csv", sep=",",header=TRUE)

samp<-paste("Rarefaction/Curves", ".pdf", sep="")
pdf(file = samp, width = 8, height = 5, bg="transparent")
par(mfrow=c(1,1), mar = c(2, 3, 0.1, 0.1))
plot(RARE1$sample, RARE1$OTUs, ylim=c(0,3000), xlim=c(0,200000), pch=19, cex=0.25, col="blue", axes=FALSE, ylab="", xlab="")
par(new=TRUE)
plot(RARE2$sample, RARE2$OTUs, ylim=c(0,3000), xlim=c(0,200000), pch=19, cex=0.25, col="green", axes=FALSE, ylab="", xlab="")
par(new=TRUE)
plot(RARE3$sample, RARE3$OTUs, ylim=c(0,3000), xlim=c(0,200000), pch=19, cex=0.25, col="black", axes=FALSE, ylab="", xlab="")
par(new=TRUE)
plot(RARE4$sample, RARE4$OTUs, ylim=c(0,3000), xlim=c(0,200000), pch=19, cex=0.25, col="red", axes=FALSE, ylab="", xlab="")
box(col="black")
grid(col="grey")
axis(2, at=, col.axis="black", las=2, cex.axis=1)
axis(1, at=, col.axis="black", las=1, cex.axis=1)
dev.off()


