
library(ggplot2)
require(grid)
library(gridExtra)
library(scales)
library(lattice)
library(Biobase)
library(reshape2)
rm(list=ls())

#######   Read in data from master folder holding all sample folders - - - - - - - - - - - - - - - -
#######   Read in data from master folder holding all sample folders - - - - - - - - - - - - - - - -
#######   Read in data from master folder holding all sample folders - - - - - - - - - - - - - - - -
#######   Read in data from master folder holding all sample folders - - - - - - - - - - - - - - - -

###############
###############
###############
###############
###############
phyto<-read.csv("CuracaoDataset.csv", header=TRUE)

##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################

populations<-c("picophyto" ,"fullSyn", "synPop1", "synPop2", "synPop3", "prochl", "hetBact")
first<-list()
second<-list()
for(I in 1:7)
{
	Ave1<-dcast(data6, Sample  + Runs + category ~ time.point, value.var=populations[I])
	Ave2<-cbind(Ave1[,1:3], "hour1u"=0, "hour2u"=log(Ave1[,5]/Ave1[,4]), "hour3u"=log(Ave1[,6]/Ave1[,5]), "hour5u"=log(Ave1[,6]/Ave1[,4]))
	Ave3<-subset(Ave2, Ave2$Runs=="1" | Ave2$Runs=="3")
	Ave4<-subset(Ave2, Ave2$Runs=="2" | Ave2$Runs=="4")
	firsty<-data.frame("pop"=populations[I], Ave3)
	secondy<-data.frame("pop"=populations[I], Ave4)
	first[[I]]<-firsty
	second[[I]]<-secondy
}
harvest1<-data.frame("diel"="first", do.call(rbind, first))
harvest2<-data.frame("diel"="second", do.call(rbind, second))
harvest3<-rbind(harvest1,harvest2)
COLOR<-c("white","#a1d99b","white","#41ab5d","white","#fd8d3c","white","#fc4e2a","white","#fd8d3c","white","#fc4e2a","white","#fd8d3c","white","#fc4e2a","white","#fd8d3c","white","#fc4e2a","white","#6baed6","white","#2171b5","white","#8c6bb1","white","#88419d")
CODES<-c("ConFP","CorFP","ConSP","CorSP","ConFFS","CorFFS","ConSFS","CorSFS","ConFS1","CorFS1","ConSS1","CorSS1","ConFS2","CorFS2","ConSS2","CorSS2","ConFS3","CorFS3","ConSS3","CorSS3","ConFPR","CorFPR","ConSPR","CorSPR","ConFH","CorFH","ConSH","CorSH")
MMgr<-subset(harvest3, category=="madracis" | category=="cont")
#MMgr<-droplevels(MMgr[!MMgr$Sample=="MM1-april28",])
PAgr<-subset(harvest3, category=="porites" | category=="cont")
#PAgr<-droplevels(PAgr[!PAgr$Sample=="PA3-april29",])
SIgr<-subset(harvest3, category=="siderastrea" | category=="cont")
HOOK<-list(MMgr, PAgr, SIgr)

############## GROWTH RATES ################
############## GROWTH RATES ################
############## GROWTH RATES ################

conts<-subset(MMgr, category=="cont")
Mcont<-aggregate(conts$hour5u, by=list(conts$pop, conts$diel), FUN=mean, na.rm=TRUE)
SDcont<-aggregate(conts$hour5u, by=list(conts$pop, conts$diel), FUN=sd, na.rm=TRUE)
Mes<-dcast(Mcont, Group.1  ~ Group.2)
SDes<-dcast(SDcont, Group.1  ~ Group.2)
write.csv(Mes, "GrowthMeans.csv")
write.csv(SDes, "GrowthSD.csv")

populations<-c("picophyto", "prochl" , "synPop1", "synPop2", "synPop3", "hetBact", "fullSyn")
ttest<-list()
for(I in 1:7)
{
	actives<-subset(conts, pop==populations[I])
	ttest[[I]]<-t.test((actives$hour5u+5) ~ as.factor(actives$diel))[[3]]
}
ttest

###################FIGURE4

samp<-paste("Figures/grazer-rates-PhytoProch", ".pdf", sep="")
pdf(file = samp, width = 2.25, height = 8, bg="transparent")
###################
par(mfrow=c(3,1), oma = c(0.1, 3.5, 0.1, 0.1), mar = c(1, 0.1, 1, 0.1))
for(I in 1:3)
{
	hok<-HOOK[[I]]
	Ave3M<-data.frame(CODES,COLOR,aggregate(hok$hour5u, by=list(hok$category, hok$diel, hok$pop), FUN=mean, na.rm=TRUE))
	Ave3SD<-data.frame(CODES,COLOR,aggregate(hok$hour5u, by=list(hok$category, hok$diel, hok$pop), FUN=sd, na.rm=TRUE))
	SpacingsA<-c("ConFP","CorFP",NA,"ConSP","CorSP")
	Ave<-droplevels(Ave3M[match(SpacingsA, Ave3M$CODES),])
	SD<-droplevels(Ave3SD[match(SpacingsA, Ave3SD$CODES),])
	NONO<-as.matrix(Ave$COLOR)
	NONO[is.na(NONO)]<-"white"
	bar<-barplot(Ave$x, width=c(1), ylim=c(-1.6, 0.75), col=NONO, axes=FALSE)
	arrows(bar,(Ave$x+(SD$x/sqrt(1))), bar,(Ave$x-(SD$x/sqrt(1))), length=0.01, angle=90, code=3,  col="black", lwd=1)
	abline(h=0, col="black", lwd=1.5)
	axis(2, at=, col.axis="black", las=2, cex.axis=0.0001)
	box(col = "black")
	#grid(col="grey", lwd=0.15, lty="dotted")
}
dev.off()

samp<-paste("Figures/grazer-rates-Proch", ".pdf", sep="")
pdf(file = samp, width = 2.25, height = 8, bg="transparent")
###################
par(mfrow=c(3,1), oma = c(0.1, 3.5, 0.1, 0.1), mar = c(1, 0.1, 1, 0.1))
for(I in 1:3)
{
	hok<-HOOK[[I]]
	Ave3M<-data.frame(CODES,COLOR,aggregate(hok$hour5u, by=list(hok$category, hok$diel, hok$pop), FUN=mean, na.rm=TRUE))
	Ave3SD<-data.frame(CODES,COLOR,aggregate(hok$hour5u, by=list(hok$category, hok$diel, hok$pop), FUN=sd, na.rm=TRUE))
	SpacingsA<-c("ConFPR","CorFPR",NA,"ConSPR","CorSPR")
	Ave<-droplevels(Ave3M[match(SpacingsA, Ave3M$CODES),])
	SD<-droplevels(Ave3SD[match(SpacingsA, Ave3SD$CODES),])
	NONO<-as.matrix(Ave$COLOR)
	NONO[is.na(NONO)]<-"white"
	bar<-barplot(Ave$x, width=c(1), ylim=c(-1.6, 0.75), col=NONO, axes=FALSE)
	arrows(bar,(Ave$x+(SD$x/sqrt(1))), bar,(Ave$x-(SD$x/sqrt(1))), length=0.01, angle=90, code=3,  col="black", lwd=1)
	abline(h=0, col="black", lwd=1.5)
	#axis(2, at=, col.axis="black", las=2, cex.axis=1.9)
	box(col = "black")
	#grid(col="grey", lwd=0.15, lty="dotted")
}
dev.off()

samp<-paste("Figures/grazer-rates-fullSyn", ".pdf", sep="")
pdf(file = samp, width = 2.25, height = 8, bg="transparent")
###################
par(mfrow=c(3,1), oma = c(0.1, 3.5, 0.1, 0.1), mar = c(1, 0.1, 1, 0.1))
for(I in 1:3)
{
	hok<-HOOK[[I]]
	Ave3M<-data.frame(CODES,COLOR,aggregate(hok$hour5u, by=list(hok$category, hok$diel, hok$pop), FUN=mean, na.rm=TRUE))
	Ave3SD<-data.frame(CODES,COLOR,aggregate(hok$hour5u, by=list(hok$category, hok$diel, hok$pop), FUN=sd, na.rm=TRUE))
	SpacingsA<-c("ConFFS","CorFFS",NA,"ConSFS","CorSFS")
	Ave<-droplevels(Ave3M[match(SpacingsA, Ave3M$CODES),])
	SD<-droplevels(Ave3SD[match(SpacingsA, Ave3SD$CODES),])
	NONO<-as.matrix(Ave$COLOR)
	NONO[is.na(NONO)]<-"white"
	par(lwd=2)
	bar<-barplot(Ave$x, width=c(1), ylim=c(-1.6, 0.75), col=NONO, axes=FALSE, lwd=0.5)
	arrows(bar,(Ave$x+(SD$x/sqrt(1))), bar,(Ave$x-(SD$x/sqrt(1))), length=0.01, angle=90, code=3,  col="black", lwd=1)
	abline(h=0, col="black", lwd=1.5)
	#axis(2, at=, col.axis="black", las=2, cex.axis=1.9)
	box(col = "black", lwd=0.75)
	#grid(col="grey", lwd=0.15, lty="dotted")
}
dev.off()

samp<-paste("Figures/grazer-rates_SYN", ".pdf", sep="")
pdf(file = samp, width = 6.75, height = 8, bg="transparent")
###################
par(mfrow=c(3,1), oma = c(0.1, 0.1, 0.1, 0.1), mar = c(1, 0.1, 1, 0.1))
for(I in 1:3)
{
	hok<-HOOK[[I]]
	Ave3M<-data.frame(CODES,COLOR,aggregate(hok$hour5u, by=list(hok$category, hok$diel, hok$pop), FUN=mean, na.rm=TRUE))
	Ave3SD<-data.frame(CODES,COLOR,aggregate(hok$hour5u, by=list(hok$category, hok$diel, hok$pop), FUN=sd, na.rm=TRUE))
	SpacingsB<-c("ConFS1","CorFS1",NA,"ConSS1","CorSS1",NA,NA,"ConFS2","CorFS2",NA,"ConSS2","CorSS2",NA,NA,"ConFS3","CorFS3",NA,"ConSS3","CorSS3")
	Ave<-droplevels(Ave3M[match(SpacingsB, Ave3M$CODES),])
	SD<-droplevels(Ave3SD[match(SpacingsB, Ave3SD$CODES),])
	NONO<-as.matrix(Ave$COLOR)
	NONO[is.na(NONO)]<-"white"
	bar<-barplot(Ave$x, width=c(1), ylim=c(-1.6, 0.75), col=NONO, axes=FALSE)
	arrows(bar,(Ave$x+(SD$x/sqrt(1))), bar,(Ave$x-(SD$x/sqrt(1))), length=0.01, angle=90, code=3,  col="black", lwd=1)
	abline(h=0, col="black", lwd=1.5)
	#axis(2, at=, col.axis="black", las=2, cex.axis=1.65)
	box(col = "black")
	#grid(col="grey", lwd=0.15, lty="dotted")
}
dev.off()


samp<-paste("Figures/grazer-rates-HET", ".pdf", sep="")
pdf(file = samp, width = 2.25, height = 8, bg="transparent")
###################
par(mfrow=c(3,1), oma = c(0.1, 0.1, 0.1, 0.1), mar = c(1, 0.1, 1, 0.1))
for(I in 1:3)
{
	hok<-HOOK[[I]]
	Ave3M<-data.frame(CODES,COLOR,aggregate(hok$hour5u, by=list(hok$category, hok$diel, hok$pop), FUN=mean, na.rm=TRUE))
	Ave3SD<-data.frame(CODES,COLOR,aggregate(hok$hour5u, by=list(hok$category, hok$diel, hok$pop), FUN=sd, na.rm=TRUE))
	SpacingsA<-c("ConFH","CorFH",NA,"ConSH","CorSH")
	Ave<-droplevels(Ave3M[match(SpacingsA, Ave3M$CODES),])
	SD<-droplevels(Ave3SD[match(SpacingsA, Ave3SD$CODES),])
	NONO<-as.matrix(Ave$COLOR)
	NONO[is.na(NONO)]<-"white"
	bar<-barplot(Ave$x, width=c(1), ylim=c(-1.6, 0.75), col=NONO , axes=FALSE)
	arrows(bar,(Ave$x+(SD$x/sqrt(1))), bar,(Ave$x-(SD$x/sqrt(1))), length=0.01, angle=90, code=3,  col="black", lwd=1)
	abline(h=0, col="black", lwd=1.5)
	#axis(2, at=, col.axis="black", las=2, cex.axis=1)
	box(col = "black")
	#grid(col="grey", lwd=0.15, lty="dotted")
}
dev.off()


##########################
##########################
##########################
##########################
library(lmerTest)
library(lme4)
library(broom)
library(multcomp)

dat2<-list()
for (D in 1:3)
{	active1<-HOOK[[D]]
	active2<-data.frame("sections"=paste(active1$diel, active1$pop, sep="-"), active1)
	leveled<-levels(active2$sections)
	numbb<-nlevels(active2$sections)
	dat1<-list()
	for (I in 1:numbb)
	{	sect<-leveled[[I]]
		active3<-subset(active2,sections==leveled[[I]])
		active4<-data.frame("Runs"=active3$Runs, "sample"=active3$Sample, "category"=active3$category, active3[,8:9])
		active5<-na.omit(melt(active4, id.vars=c("Runs", "sample", "category"), variable.name="time", value.name="change"))
		active5$sample<-droplevels(active5$sample)
		active5$time<-droplevels(active5$time)
		fit1 <- lmer(change ~ category*time + (1|Runs/time) + (1|Runs/sample), data=active5, REML=FALSE)
		tmSer<-anova(fit1)[[6]][1]
		spec<-active3[5,6]
		ttest<-t.test((active3$hour5u+5) ~ active3$category)[[3]]
		Wilcoxon<-wilcox.test((active3$hour5u+5) ~ factor(active3$category), p.adjust="none", paired=FALSE)[[3]]
		shap<-shapiro.test(active3$hour5u)[[2]]
		dat1[[I]]<-data.frame(spec, sect, round(shap,4), round(ttest,4), round(Wilcoxon,4), round(tmSer,4))
	}
	dat2[[D]]<-data.frame(do.call(rbind, dat1))
}
dat2


sets<-rbind(dat2[[1]], dat2[[2]], dat2[[3]])

ss<-subset(sets, sect=="first-picophyto")
data.frame(ss[,1:2], "Padjust"=round(p.adjust(ss[,6], method="fdr", n=length(ss[,6])), 4))

ss<-subset(sets, sect=="second-picophyto")
data.frame(ss[,1:2], "Padjust"=round(p.adjust(ss[,6], method="fdr", n=length(ss[,6])), 4))

ss<-subset(sets, sect=="first-prochl")
data.frame(ss[,1:2], "Padjust"=round(p.adjust(ss[,6], method="bonferroni", n=length(ss[,6])), 4))

ss<-subset(sets, sect=="second-prochl")
data.frame(ss[,1:2], "Padjust"=round(p.adjust(ss[,6], method="fdr", n=length(ss[,6])), 4))

ss<-subset(sets, sect=="first-fullSyn")
data.frame(ss[,1:2], "Padjust"=round(p.adjust(ss[,6], method="fdr", n=length(ss[,6])), 4))

ss<-subset(sets, sect=="second-fullSyn")
data.frame(ss[,1:2], "Padjust"=round(p.adjust(ss[,6], method="fdr", n=length(ss[,6])), 4))

ss<-subset(sets, sect=="first-synPop1")
data.frame(ss[,1:2], "Padjust"=round(p.adjust(ss[,6], method="fdr", n=length(ss[,6])), 4))

ss<-subset(sets, sect=="second-synPop1")
data.frame(ss[,1:2], "Padjust"=round(p.adjust(ss[,6], method="fdr", n=length(ss[,6])), 4))

ss<-subset(sets, sect=="first-synPop2")
data.frame(ss[,1:2], "Padjust"=round(p.adjust(ss[,6], method="fdr", n=length(ss[,6])), 4))

ss<-subset(sets, sect=="second-synPop2")
data.frame(ss[,1:2], "Padjust"=round(p.adjust(ss[,6], method="fdr", n=length(ss[,6])), 4))

ss<-subset(sets, sect=="first-synPop3")
data.frame(ss[,1:2], "Padjust"=round(p.adjust(ss[,6], method="fdr", n=length(ss[,6])), 4))

ss<-subset(sets, sect=="second-synPop3")
data.frame(ss[,1:2], "Padjust"=round(p.adjust(ss[,6], method="fdr", n=length(ss[,6])), 4))

ss<-subset(sets, sect=="first-hetBact")
data.frame(ss[,1:2], "Padjust"=round(p.adjust(ss[,6], method="fdr", n=length(ss[,6])), 4))

ss<-subset(sets, sect=="second-hetBact")
data.frame(ss[,1:2], "Padjust"=round(p.adjust(ss[,6], method="none", n=length(ss[,6])), 4))

write.csv(dat2, "stats-figure4.csv")


###########
###########
###########
###########
###########
###########
picoCarb<-530
picoNitro<-1
picoPhos<-1
ProCarb<-39
ProNitro<-9.6
ProPhos<-0.34
SynCarb<-82
SynNitro<-39.8
SynPhos<-0.81
Nitro<-c(SynNitro, ProNitro, picoNitro)
Carbo<-c(SynCarb, ProCarb, picoCarb)
Phos<-c(SynPhos, ProPhos, picoPhos)
populations<-c("fullSyn", "prochl", "picophyto")
means<-list()
SDs<-list()
datsta<-list()
st.err <- function(x) {sd(x)/sqrt(length(x))} 
for(I in 1:3)
{
	popu<-populations[I]

	Ave1<-dcast(data6, Sample  + Runs + category ~ time.point, value.var=popu)
	R1<-subset(Ave1, Runs=="1")
	ContR1<-subset(aggregate(R1[,4:6], by=list(R1[,3]), FUN=mean, na.rm=TRUE), Group.1=="cont")
	AveR1<-cbind(R1[,1:3], "0"= (ContR1[,2]-R1[,4]), "1"=(ContR1[,3]-R1[,5]), "5N"=((ContR1[,4]-R1[,6])*400*Nitro[I]), "5C"=((ContR1[,4]-R1[,6])*400*Carbo[I]), "5P"=((ContR1[,4]-R1[,6])*400*Phos[I]))
	#R2<-subset(Ave1, Runs=="2")
	#ContR2<-subset(aggregate(R2[,4:6], by=list(R2[,3]), FUN=mean, na.rm=TRUE), Group.1=="cont")
	#AveR2<-cbind(R2[,1:3], "0"= (ContR2[,2]-R2[,4]), "1"=(ContR2[,3]-R2[,5]), "5N"=((ContR2[,4]-R2[,6])*400*Nitro[I]), "5C"=((ContR2[,4]-R2[,6])*400*Carbo[I]), "5P"=((ContR2[,4]-R2[,6])*400*Phos[I]))
	R3<-subset(Ave1, Runs=="3")
	ContR3<-subset(aggregate(R3[,4:6], by=list(R3[,3]), FUN=mean, na.rm=TRUE), Group.1=="cont")
	AveR3<-cbind(R3[,1:3], "0"= (ContR3[,2]-R3[,4]), "1"=(ContR3[,3]-R3[,5]), "5N"=((ContR3[,4]-R3[,6])*400*Nitro[I]), "5C"=((ContR3[,4]-R3[,6])*400*Carbo[I]), "5P"=((ContR3[,4]-R3[,6])*400*Phos[I]))
	#R4<-subset(Ave1, Runs=="4")
	#ContR4<-subset(aggregate(R4[,4:6], by=list(R4[,3]), FUN=mean, na.rm=TRUE), Group.1=="cont")
	#AveR4<-cbind(R4[,1:3], "0"= (ContR4[,2]-R4[,4]), "1"=(ContR4[,3]-R4[,5]), "5N"=((ContR4[,4]-R4[,6])*400*Nitro[I]), "5C"=((ContR4[,4]-R4[,6])*400*Carbo[I]), "5P"=((ContR4[,4]-R4[,6])*400*Phos[I]))
	
	Ave5<-rbind(AveR1, AveR3)
	Ave5
	
	#Ave5<-Ave5[-32,]
	#Ave5<-Ave5[-16,]
	#Ave5<-Ave5[-14,]
	#Ave5<-Ave5[-9:-11,]
	datsta[[I]]<-Ave5
	means[[I]]<-aggregate(Ave5[,6:8]/5, by=list(Ave5$category), FUN=mean, na.rm=TRUE)
	SDs[[I]]<-aggregate(Ave5[,6:8]/5, by=list(Ave5$category), FUN=sd, na.rm=TRUE)

}

ddd<-data.frame(means[[2]][2:4,3], means[[2]][2:4,2], means[[2]][2:4,4])
colnames(ddd)<-c("carbon", "nitrogen", "phosphate")
rownames(ddd)<-c("madracis", "porites", "siderastrea")
dddd<-ddd/1000000
t(dddd)

ddd<-data.frame(SDs[[2]][2:4,3], SDs[[2]][2:4,2], SDs[[2]][2:4,4])
colnames(ddd)<-c("carbon", "nitrogen", "phosphate")
rownames(ddd)<-c("madracis", "porites", "siderastrea")
dddd<-ddd/1000000
t(dddd)

ddd<-data.frame(means[[1]][2:4,3], means[[1]][2:4,2], means[[1]][2:4,4])
colnames(ddd)<-c("carbon", "nitrogen", "phosphate")
rownames(ddd)<-c("madracis", "porites", "siderastrea")
dddd<-ddd/1000000
t(dddd)

ddd<-data.frame(SDs[[1]][2:4,3], SDs[[1]][2:4,2], SDs[[1]][2:4,4])
colnames(ddd)<-c("carbon", "nitrogen", "phosphate")
rownames(ddd)<-c("madracis", "porites", "siderastrea")
dddd<-ddd/1000000
t(dddd)

ddd<-data.frame(means[[3]][2:4,3], means[[3]][2:4,2], means[[3]][2:4,4])
colnames(ddd)<-c("carbon", "nitrogen", "phosphate")
rownames(ddd)<-c("madracis", "porites", "siderastrea")
dddd<-ddd/1000000
t(dddd)

ddd<-data.frame(SDs[[3]][2:4,3], SDs[[3]][2:4,2], SDs[[3]][2:4,4])
colnames(ddd)<-c("carbon", "nitrogen", "phosphate")
rownames(ddd)<-c("madracis", "porites", "siderastrea")
dddd<-ddd/1000000
t(dddd)

one<-data.frame(datsta[1])
two<-data.frame(datsta[2])
three<-data.frame(datsta[3])

combi<-cbind(one[,1:3],one[,7], two[,7], three[,7])
com<-data.frame(one[,1:3],data.frame(rowSums(combi[,4:6], na.rm=TRUE)))
aggregate((com[,4]/5)/1000000, by=list(com$category), FUN=mean, na.rm=TRUE)
aggregate((com[,4]/5)/1000000, by=list(com$category), FUN=sd, na.rm=TRUE)

############
############

visualizer<-function(data6)
{
samp<-paste("Figures/Controls", ".pdf", sep="")
pdf(file = samp, width = 2.5, height = 6, bg="transparent")
	
populations<-c("fullSyn", "synPop1", "synPop2", "synPop3", "prochl", "picophyto", "hetBact")
par(mfrow=c(6,2), oma = c(2.0, 5, 0.1, 0.1), mar = c(0.2, 0.1, 0.2, 0.1))

R1<-subset(data6, Runs=="1")
R2<-subset(data6, Runs=="2")
R3<-subset(data6, Runs=="3")
R4<-subset(data6, Runs=="4")

popu<-populations[6]
Ave1<-dcast(R1, Sample  + Runs + category ~ time.point, value.var=popu)
plot(c(0,1,5), Ave1[1,4:6], pch=19, cex=1, col="grey", ylim=c(0, 6000), axes=FALSE)
lines(c(0,1,5), Ave1[1,4:6], col="grey", lwd=1.5)
par(new="TRUE")
plot(c(-0.1,0.9,4.9), Ave1[2,4:6], pch=19, cex=1, col="grey", ylim=c(0, 6000), axes=FALSE)
lines(c(-0.1,0.9,4.9), Ave1[2,4:6], col="grey", lwd=1.5)
Ave2<-dcast(R3, Sample  + Runs + category ~ time.point, value.var=popu)
par(new="TRUE")
plot(c(0.1,1.1,5.1), Ave2[1,4:6], pch=19, cex=1, col="grey", ylim=c(0, 6000), axes=FALSE)
lines(c(0.1,1.1,5.1), Ave2[1,4:6], col="grey", lwd=1.5)
par(new="TRUE")
plot(c(0.1,1.1,5.1), Ave2[2,4:6], pch=19, cex=1, col="grey", ylim=c(0, 6000), axes=FALSE)
lines(c(0.1,1.1,5.1), Ave2[2,4:6], col="grey", lwd=1.5)
box(col = "black")
axis(2, at=seq(0,5e+03, by=2e+03), col.axis="black", las=2, cex.axis=1.2)


Ave3<-dcast(R2, Sample  + Runs + category ~ time.point, value.var=popu)
plot(c(0,1,5), Ave3[1,4:6], pch=19, cex=1, col="black", ylim=c(0, 6000), axes=FALSE)
lines(c(0,1,5), Ave3[1,4:6], col="black", lwd=1.5)
par(new="TRUE")
plot(c(-0.1,0.9,4.9), Ave3[2,4:6], pch=19, cex=1, col="black", ylim=c(0, 6000), axes=FALSE)
lines(c(-0.1,0.9,4.9), Ave3[2,4:6], col="black", lwd=1.5)
Ave4<-dcast(R4, Sample  + Runs + category ~ time.point, value.var=popu)
par(new="TRUE")
plot(c(0.1,1.1,5.1), Ave4[1,4:6], pch=19, cex=1, col="black", ylim=c(0, 6000), axes=FALSE)
lines(c(0.1,1.1,5.1), Ave4[1,4:6], col="black", lwd=1.5)
par(new="TRUE")
plot(c(0.1,1.1,5.1), Ave4[2,4:6], pch=19, cex=1, col="black", ylim=c(0, 6000), axes=FALSE)
lines(c(0.1,1.1,5.1), Ave4[2,4:6], col="black", lwd=1.5)
box(col = "black")


popu<-populations[5]
Ave1<-dcast(R1, Sample  + Runs + category ~ time.point, value.var=popu)
plot(c(0,1,5), Ave1[1,4:6], pch=19, cex=1, col="grey", ylim=c(0, 200000), axes=FALSE)
lines(c(0,1,5), Ave1[1,4:6], col="grey", lwd=1.5)
par(new="TRUE")
plot(c(-0.1,0.9,4.9), Ave1[2,4:6], pch=19, cex=1, col="grey", ylim=c(0, 200000), axes=FALSE)
lines(c(-0.1,0.9,4.9), Ave1[2,4:6], col="grey", lwd=1.5)
Ave2<-dcast(R3, Sample  + Runs + category ~ time.point, value.var=popu)
par(new="TRUE")
plot(c(0.1,1.1,5.1), Ave2[1,4:6], pch=19, cex=1, col="grey", ylim=c(0, 200000), axes=FALSE)
lines(c(0.1,1.1,5.1), Ave2[1,4:6], col="grey", lwd=1.5)
par(new="TRUE")
plot(c(0.1,1.1,5.1), Ave2[2,4:6], pch=19, cex=1, col="grey", ylim=c(0, 200000), axes=FALSE)
lines(c(0.1,1.1,5.1), Ave2[2,4:6], col="grey", lwd=1.5)
box(col = "black")
axis(2, at=seq(0,15e+04, by=5e+04), col.axis="black", las=2, cex.axis=1.2)


Ave3<-dcast(R2, Sample  + Runs + category ~ time.point, value.var=popu)
plot(c(0,1,5), Ave3[1,4:6], pch=19, cex=1, col="black", ylim=c(0, 200000), axes=FALSE)
lines(c(0,1,5), Ave3[1,4:6], col="black", lwd=1.5)
par(new="TRUE")
plot(c(-0.1,0.9,4.9), Ave3[2,4:6], pch=19, cex=1, col="black", ylim=c(0, 200000), axes=FALSE)
lines(c(-0.1,0.9,4.9), Ave3[2,4:6], col="black", lwd=1.5)
Ave4<-dcast(R4, Sample  + Runs + category ~ time.point, value.var=popu)
par(new="TRUE")
plot(c(0.1,1.1,5.1), Ave4[1,4:6], pch=19, cex=1, col="black", ylim=c(0, 200000), axes=FALSE)
lines(c(0.1,1.1,5.1), Ave4[1,4:6], col="black", lwd=1.5)
par(new="TRUE")
plot(c(0.1,1.1,5.1), Ave4[2,4:6], pch=19, cex=1, col="black", ylim=c(0, 200000), axes=FALSE)
lines(c(0.1,1.1,5.1), Ave4[2,4:6], col="black", lwd=1.5)
box(col = "black")


popu<-populations[2]
Ave1<-dcast(R1, Sample  + Runs + category ~ time.point, value.var=popu)
plot(c(0,1,5), Ave1[1,4:6], pch=19, cex=1, col="grey", ylim=c(0, 15000), axes=FALSE)
lines(c(0,1,5), Ave1[1,4:6], col="grey", lwd=1.5)
par(new="TRUE")
plot(c(-0.1,0.9,4.9), Ave1[2,4:6], pch=19, cex=1, col="grey", ylim=c(0, 15000), axes=FALSE)
lines(c(-0.1,0.9,4.9), Ave1[2,4:6], col="grey", lwd=1.5)
Ave2<-dcast(R3, Sample  + Runs + category ~ time.point, value.var=popu)
par(new="TRUE")
plot(c(0.1,1.1,5.1), Ave2[1,4:6], pch=19, cex=1, col="grey", ylim=c(0, 15000), axes=FALSE)
lines(c(0.1,1.1,5.1), Ave2[1,4:6], col="grey", lwd=1.5)
par(new="TRUE")
plot(c(0.1,1.1,5.1), Ave2[2,4:6], pch=19, cex=1, col="grey", ylim=c(0, 15000), axes=FALSE)
lines(c(0.1,1.1,5.1), Ave2[2,4:6], col="grey", lwd=1.5)
box(col = "black")
axis(2, at=seq(0,15e+03, by=5e+03), col.axis="black", las=2, cex.axis=1.2)


Ave3<-dcast(R2, Sample  + Runs + category ~ time.point, value.var=popu)
plot(c(0,1,5), Ave3[1,4:6], pch=19, cex=1, col="black", ylim=c(0, 15000), axes=FALSE)
lines(c(0,1,5), Ave3[1,4:6], col="black", lwd=1.5)
par(new="TRUE")
plot(c(-0.1,0.9,4.9), Ave3[2,4:6], pch=19, cex=1, col="black", ylim=c(0, 15000), axes=FALSE)
lines(c(-0.1,0.9,4.9), Ave3[2,4:6], col="black", lwd=1.5)
Ave4<-dcast(R4, Sample  + Runs + category ~ time.point, value.var=popu)
par(new="TRUE")
plot(c(0.1,1.1,5.1), Ave4[1,4:6], pch=19, cex=1, col="black", ylim=c(0, 15000), axes=FALSE)
lines(c(0.1,1.1,5.1), Ave4[1,4:6], col="black", lwd=1.5)
par(new="TRUE")
plot(c(0.1,1.1,5.1), Ave4[2,4:6], pch=19, cex=1, col="black", ylim=c(0, 15000), axes=FALSE)
lines(c(0.1,1.1,5.1), Ave4[2,4:6], col="black", lwd=1.5)
box(col = "black")


popu<-populations[3]
Ave1<-dcast(R1, Sample  + Runs + category ~ time.point, value.var=popu)
plot(c(0,1,5), Ave1[1,4:6], pch=19, cex=1, col="grey", ylim=c(0, 15000), axes=FALSE)
lines(c(0,1,5), Ave1[1,4:6], col="grey", lwd=1.5)
par(new="TRUE")
plot(c(-0.1,0.9,4.9), Ave1[2,4:6], pch=19, cex=1, col="grey", ylim=c(0, 15000), axes=FALSE)
lines(c(-0.1,0.9,4.9), Ave1[2,4:6], col="grey", lwd=1.5)
Ave2<-dcast(R3, Sample  + Runs + category ~ time.point, value.var=popu)
par(new="TRUE")
plot(c(0.1,1.1,5.1), Ave2[1,4:6], pch=19, cex=1, col="grey", ylim=c(0, 15000), axes=FALSE)
lines(c(0.1,1.1,5.1), Ave2[1,4:6], col="grey", lwd=1.5)
par(new="TRUE")
plot(c(0.1,1.1,5.1), Ave2[2,4:6], pch=19, cex=1, col="grey", ylim=c(0, 15000), axes=FALSE)
lines(c(0.1,1.1,5.1), Ave2[2,4:6], col="grey", lwd=1.5)
box(col = "black")
axis(2, at=seq(0,15e+03, by=5e+03), col.axis="black", las=2, cex.axis=1.2)


Ave3<-dcast(R2, Sample  + Runs + category ~ time.point, value.var=popu)
plot(c(0,1,5), Ave3[1,4:6], pch=19, cex=1, col="black", ylim=c(0, 15000), axes=FALSE)
lines(c(0,1,5), Ave3[1,4:6], col="black", lwd=1.5)
par(new="TRUE")
plot(c(-0.1,0.9,4.9), Ave3[2,4:6], pch=19, cex=1, col="black", ylim=c(0, 15000), axes=FALSE)
lines(c(-0.1,0.9,4.9), Ave3[2,4:6], col="black", lwd=1.5)
Ave4<-dcast(R4, Sample  + Runs + category ~ time.point, value.var=popu)
par(new="TRUE")
plot(c(0.1,1.1,5.1), Ave4[1,4:6], pch=19, cex=1, col="black", ylim=c(0, 15000), axes=FALSE)
lines(c(0.1,1.1,5.1), Ave4[1,4:6], col="black", lwd=1.5)
par(new="TRUE")
plot(c(0.1,1.1,5.1), Ave4[2,4:6], pch=19, cex=1, col="black", ylim=c(0, 15000), axes=FALSE)
lines(c(0.1,1.1,5.1), Ave4[2,4:6], col="black", lwd=1.5)
box(col = "black")


popu<-populations[4]
Ave1<-dcast(R1, Sample  + Runs + category ~ time.point, value.var=popu)
plot(c(0,1,5), Ave1[1,4:6], pch=19, cex=1, col="grey", ylim=c(0, 15000), axes=FALSE)
lines(c(0,1,5), Ave1[1,4:6], col="grey", lwd=1.5)
par(new="TRUE")
plot(c(-0.1,0.9,4.9), Ave1[2,4:6], pch=19, cex=1, col="grey", ylim=c(0, 15000), axes=FALSE)
lines(c(-0.1,0.9,4.9), Ave1[2,4:6], col="grey", lwd=1.5)
Ave2<-dcast(R3, Sample  + Runs + category ~ time.point, value.var=popu)
par(new="TRUE")
plot(c(0.1,1.1,5.1), Ave2[1,4:6], pch=19, cex=1, col="grey", ylim=c(0, 15000), axes=FALSE)
lines(c(0.1,1.1,5.1), Ave2[1,4:6], col="grey", lwd=1.5)
par(new="TRUE")
plot(c(0.1,1.1,5.1), Ave2[2,4:6], pch=19, cex=1, col="grey", ylim=c(50, 15000), axes=FALSE)
lines(c(0.1,1.1,5.1), Ave2[2,4:6], col="grey", lwd=1.5)
box(col = "black")
axis(2, at=seq(0,15e+03, by=5e+03), col.axis="black", las=2, cex.axis=1.2)


Ave3<-dcast(R2, Sample  + Runs + category ~ time.point, value.var=popu)
plot(c(0,1,5), Ave3[1,4:6], pch=19, cex=1, col="black", ylim=c(0, 15000), axes=FALSE)
lines(c(0,1,5), Ave3[1,4:6], col="black", lwd=1.5)
par(new="TRUE")
plot(c(-0.1,0.9,4.9), Ave3[2,4:6], pch=19, cex=1, col="black", ylim=c(0, 15000), axes=FALSE)
lines(c(-0.1,0.9,4.9), Ave3[2,4:6], col="black", lwd=1.5)
Ave4<-dcast(R4, Sample  + Runs + category ~ time.point, value.var=popu)
par(new="TRUE")
plot(c(0.1,1.1,5.1), Ave4[1,4:6], pch=19, cex=1, col="black", ylim=c(0, 15000), axes=FALSE)
lines(c(0.1,1.1,5.1), Ave4[1,4:6], col="black", lwd=1.5)
par(new="TRUE")
plot(c(0.1,1.1,5.1), Ave4[2,4:6], pch=19, cex=1, col="black", ylim=c(0, 15000), axes=FALSE)
lines(c(0.1,1.1,5.1), Ave4[2,4:6], col="black", lwd=1.5)
box(col = "black")


popu<-populations[7]
Ave1<-dcast(R1, Sample  + Runs + category ~ time.point, value.var=popu)
plot(c(0,1,5), Ave1[1,4:6], pch=19, cex=1, col="grey", ylim=c(0, 800000), axes=FALSE)
lines(c(0,1,5), Ave1[1,4:6], col="grey", lwd=1.5)
par(new="TRUE")
plot(c(-0.1,0.9,4.9), Ave1[2,4:6], pch=19, cex=1, col="grey", ylim=c(0, 800000), axes=FALSE)
lines(c(-0.1,0.9,4.9), Ave1[2,4:6], col="grey", lwd=1.5)
Ave2<-dcast(R3, Sample  + Runs + category ~ time.point, value.var=popu)
par(new="TRUE")
plot(c(0.1,1.1,5.1), Ave2[1,4:6], pch=19, cex=1, col="grey", ylim=c(0, 800000), axes=FALSE)
lines(c(0.1,1.1,5.1), Ave2[1,4:6], col="grey", lwd=1.5)
par(new="TRUE")
plot(c(0.1,1.1,5.1), Ave2[2,4:6], pch=19, cex=1, col="grey", ylim=c(0, 800000), axes=FALSE)
lines(c(0.1,1.1,5.1), Ave2[2,4:6], col="grey", lwd=1.5)
box(col = "black")
axis(1, at=, col.axis="black", las=1, cex.axis=1.2)
axis(2, at=seq(0,8e+05, by=2e+05), col.axis="black", las=1, cex.axis=1.2)

Ave3<-dcast(R2, Sample  + Runs + category ~ time.point, value.var=popu)
plot(c(0,1,5), Ave3[1,4:6], pch=19, cex=1, col="black", ylim=c(0, 800000), axes=FALSE)
lines(c(0,1,5), Ave3[1,4:6], col="black", lwd=1.5)
par(new="TRUE")
plot(c(-0.1,0.9,4.9), Ave3[2,4:6], pch=19, cex=1, col="black", ylim=c(0, 800000), axes=FALSE)
lines(c(-0.1,0.9,4.9), Ave3[2,4:6], col="black", lwd=1.5)
Ave4<-dcast(R4, Sample  + Runs + category ~ time.point, value.var=popu)
par(new="TRUE")
plot(c(0.1,1.1,5.1), Ave4[1,4:6], pch=19, cex=1, col="black", ylim=c(0, 800000), axes=FALSE)
lines(c(0.1,1.1,5.1), Ave4[1,4:6], col="black", lwd=1.5)
par(new="TRUE")
plot(c(0.1,1.1,5.1), Ave4[2,4:6], pch=19, cex=1, col="black", ylim=c(0, 800000), axes=FALSE)
lines(c(0.1,1.1,5.1), Ave4[2,4:6], col="black", lwd=1.5)
box(col = "black")
axis(1, at=, col.axis="black", las=1, cex.axis=1.2)

dev.off()

}

	
populations<-c("fullSyn", "synPop1", "synPop2", "synPop3", "prochl", "picophyto", "hetBact")

R1<-subset(data6, Runs=="1")
R2<-subset(data6, Runs=="2")
R3<-subset(data6, Runs=="3")
R4<-subset(data6, Runs=="4")

popu<-populations[6]
Ave1<-dcast(R1, Sample  + Runs + category ~ time.point, value.var=popu)
Ave2<-dcast(R3, Sample  + Runs + category ~ time.point, value.var=popu)
Avve<-rbind(Ave1, Ave2)
AVVVe1<-data.frame(subset(Avve, category=="cont"), popu)

popu<-populations[5]
Ave1<-dcast(R1, Sample  + Runs + category ~ time.point, value.var=popu)
Ave2<-dcast(R3, Sample  + Runs + category ~ time.point, value.var=popu)
Avve<-rbind(Ave1, Ave2)
AVVVe2<-data.frame(subset(Avve, category=="cont"), popu)

popu<-populations[1]
Ave1<-dcast(R1, Sample  + Runs + category ~ time.point, value.var=popu)
Ave2<-dcast(R3, Sample  + Runs + category ~ time.point, value.var=popu)
Avve<-rbind(Ave1, Ave2)
AVVVe3<-data.frame(subset(Avve, category=="cont"), popu)

popu<-populations[2]
Ave1<-dcast(R1, Sample  + Runs + category ~ time.point, value.var=popu)
Ave2<-dcast(R3, Sample  + Runs + category ~ time.point, value.var=popu)
Avve<-rbind(Ave1, Ave2)
AVVVe4<-data.frame(subset(Avve, category=="cont"), popu)

popu<-populations[3]
Ave1<-dcast(R1, Sample  + Runs + category ~ time.point, value.var=popu)
Ave2<-dcast(R3, Sample  + Runs + category ~ time.point, value.var=popu)
Avve<-rbind(Ave1, Ave2)
AVVVe5<-data.frame(subset(Avve, category=="cont"), popu)

popu<-populations[4]
Ave1<-dcast(R1, Sample  + Runs + category ~ time.point, value.var=popu)
Ave2<-dcast(R3, Sample  + Runs + category ~ time.point, value.var=popu)
Avve<-rbind(Ave1, Ave2)
AVVVe6<-data.frame(subset(Avve, category=="cont"), popu)

popu<-populations[7]
Ave1<-dcast(R1, Sample  + Runs + category ~ time.point, value.var=popu)
Ave2<-dcast(R3, Sample  + Runs + category ~ time.point, value.var=popu)
Avve<-rbind(Ave1, Ave2)
AVVVe7<-data.frame(subset(Avve, category=="cont"), popu)

postsunset<-rbind(AVVVe1, AVVVe2, AVVVe3, AVVVe4, AVVVe5, AVVVe6, AVVVe7)
psMean1<-aggregate(postsunset$X0, by=list(postsunset$popu), FUN=mean, na.rm=TRUE)
psSD1<-aggregate(postsunset$X0, by=list(postsunset$popu), FUN=sd, na.rm=TRUE)

fgCar<-rbind(530,39,82,82,82,82,20)
psCarMean1<-data.frame(psMean1$Group.1, "x"=(psMean1$x*fgCar[,1])/1000000)
psCarSD1<-data.frame(psSD1$Group.1, "x"=(psSD1$x*fgCar[,1])/1000000)

popu<-populations[6]
Ave1<-dcast(R2, Sample  + Runs + category ~ time.point, value.var=popu)
Ave2<-dcast(R4, Sample  + Runs + category ~ time.point, value.var=popu)
Avve<-rbind(Ave1, Ave2)
AVVVe1<-data.frame(subset(Avve, category=="cont"), popu)

popu<-populations[5]
Ave1<-dcast(R2, Sample  + Runs + category ~ time.point, value.var=popu)
Ave2<-dcast(R4, Sample  + Runs + category ~ time.point, value.var=popu)
Avve<-rbind(Ave1, Ave2)
AVVVe2<-data.frame(subset(Avve, category=="cont"), popu)

popu<-populations[1]
Ave1<-dcast(R2, Sample  + Runs + category ~ time.point, value.var=popu)
Ave2<-dcast(R4, Sample  + Runs + category ~ time.point, value.var=popu)
Avve<-rbind(Ave1, Ave2)
AVVVe3<-data.frame(subset(Avve, category=="cont"), popu)

popu<-populations[2]
Ave1<-dcast(R2, Sample  + Runs + category ~ time.point, value.var=popu)
Ave2<-dcast(R4, Sample  + Runs + category ~ time.point, value.var=popu)
Avve<-rbind(Ave1, Ave2)
AVVVe4<-data.frame(subset(Avve, category=="cont"), popu)

popu<-populations[3]
Ave1<-dcast(R2, Sample  + Runs + category ~ time.point, value.var=popu)
Ave2<-dcast(R4, Sample  + Runs + category ~ time.point, value.var=popu)
Avve<-rbind(Ave1, Ave2)
AVVVe5<-data.frame(subset(Avve, category=="cont"), popu)

popu<-populations[4]
Ave1<-dcast(R2, Sample  + Runs + category ~ time.point, value.var=popu)
Ave2<-dcast(R4, Sample  + Runs + category ~ time.point, value.var=popu)
Avve<-rbind(Ave1, Ave2)
AVVVe6<-data.frame(subset(Avve, category=="cont"), popu)

popu<-populations[7]
Ave1<-dcast(R2, Sample  + Runs + category ~ time.point, value.var=popu)
Ave2<-dcast(R4, Sample  + Runs + category ~ time.point, value.var=popu)
Avve<-rbind(Ave1, Ave2)
AVVVe7<-data.frame(subset(Avve, category=="cont"), popu)

presunrise<-rbind(AVVVe1, AVVVe2, AVVVe3, AVVVe4, AVVVe5, AVVVe6, AVVVe7)
#postsunset<-rbind(AVVVe1, AVVVe2, AVVVe3, AVVVe4, AVVVe5, AVVVe6, AVVVe7)
psMean2<-aggregate(presunrise$X0, by=list(presunrise$popu), FUN=mean, na.rm=TRUE)
psSD2<-aggregate(presunrise$X0, by=list(presunrise$popu), FUN=sd, na.rm=TRUE)

fgCar<-rbind(530,39,82,82,82,82,20)
psCarMean2<-data.frame(psMean2$Group.1, "x"=(psMean2$x*fgCar[,1])/1000000)
psCarSD2<-data.frame(psSD2$Group.1, "x"=(psSD2$x*fgCar[,1])/1000000)

samp<-paste("Figures/ControlsBARS", ".pdf", sep="")
pdf(file = samp, width = 8, height = 4, bg="transparent")

par(mfrow=c(1,1), oma = c(0.1, 0.1, 0.1, 0.1), mar = c(0.75, 0.1, 0.75, 2))
psMean<-t(cbind(psMean1$x, psMean2$x))
psSD<-t(cbind(psSD1$x, psSD2$x))
bar<-barplot(psMean, beside=T, width=c(0.78),axes=FALSE, ylim=c(1000, 1000000), log="y", col=c("#f7f7f7","#f7f7f7","#d9d9d9","#d9d9d9","#969696","#969696","#fee8c8","#fee8c8","#fdbb84","#fdbb84","#fc8d59","#fc8d59","#525252","#525252"))
arrows(bar, psMean+psSD, bar, psMean, lwd=1.0, angle=90, code=3, length=0.025)
axis(2, col.axis="black", las=2, cex.axis=0.8, at=c(1000,10000,100000,1000000))
box(col="black")
grid(col="grey")

dev.off()




samp<-paste("Figures/ControlsCarbonBARS", ".pdf", sep="")
pdf(file = samp, width = 10, height = 2.5, bg="transparent")

par(mfrow=c(1,1), oma = c(0.1, 3.2, 0.1, 0.1), mar = c(0.75, 0.1, 0.75, 2))
psMean<-t(cbind(psCarMean1$x, psCarMean2$x))
psSD<-t(cbind(psCarSD1$x, psCarSD2$x))
bar<-barplot(psMean, beside=T, width=c(0.78),axes=FALSE, ylim=c(0.1, 10), log="y", col=c("#f7f7f7","#f7f7f7","#d9d9d9","#d9d9d9","#969696","#969696","#fee8c8","#fee8c8","#fdbb84","#fdbb84","#fc8d59","#fc8d59","#525252","#525252"))
arrows(bar, psMean+psSD, bar, psMean, lwd=1.0, angle=90, code=3, length=0.025)
axis(2, col.axis="black", las=2, cex.axis=1.25, at=c(0.1,1,10))
box(col="black")
grid(col="grey")

dev.off()




bbb<-c(0.4444,0.4444,0.4444,0.4444,0.5516,0.5516,0.5516,0.5516,0.2991,0.2991,0.2991,0.2991,0.5113,0.5113,0.5113,0.5113,0.2665,0.2665,0.2665,0.2665,0.1011,0.1011,0.1011,0.1011,-0.0286,-0.0286,-0.0286,-0.0286)

fff<-c(530,530,530,530,39,39,39,39,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,20,20,20,20)


prod<-data.frame(postsunset, bbb, fff)

Finprod<-data.frame(prod, "newcell"=((exp(prod$bbb)*prod$X0)-prod$X0), "ug"=((((exp(prod$bbb)*prod$X0)-prod$X0)*prod$fff)/1000))

NewCellMean<-aggregate(Finprod$newcell, by=list(Finprod$popu), FUN=mean, na.rm=TRUE)
NewCellSD<-aggregate(Finprod$newcell, by=list(Finprod$popu), FUN=sd, na.rm=TRUE)

NewCarbonMean<-aggregate(Finprod$ug, by=list(Finprod$popu), FUN=mean, na.rm=TRUE)
NewCarbonSD<-aggregate(Finprod$ug, by=list(Finprod$popu), FUN=sd, na.rm=TRUE)




post<-data.frame(postsunset, "time"="post")
pre<-data.frame(presunrise, "time"="pre")
combined<-rbind(post, pre)
picophyto<-subset(combined, popu=="picophyto")
prochlo<-subset(combined, popu=="prochl")
fullSyn<-subset(combined, popu=="fullSyn")
synPop1<-subset(combined, popu=="synPop1")
synPop2<-subset(combined, popu=="synPop2")
synPop3<-subset(combined, popu=="synPop3")
hetBact<-subset(combined, popu=="hetBact")


t.test(picophyto$X0 ~ picophyto$time)
shapiro.test(picophyto$X0)
t.test(prochlo$X0 ~ prochlo$time)
shapiro.test(prochlo$X0)
t.test(fullSyn$X0 ~ fullSyn$time)
shapiro.test(fullSyn$X0)
t.test(synPop1$X0 ~ synPop1$time)
shapiro.test(synPop1$X0)
t.test(synPop2$X0 ~ synPop2$time)
shapiro.test(synPop2$X0)
t.test(synPop3$X0 ~ synPop3$time)
shapiro.test(synPop3$X0)
t.test(hetBact$X0 ~ hetBact$time)
shapiro.test(hetBact$X0)


populations<-c("picophyto", "prochl" , "synPop1", "synPop2", "synPop3", "hetBact", "fullSyn")
means<-list()
SDs<-list()

for(I in 1:7)
{
	popu<-populations[I]
	Ave1<-dcast(data6, Sample  + Runs + category ~ time.point, value.var=popu)
	R1<-subset(Ave1, Runs=="1" | Runs=="3")
	#R1<-subset(Ave1, Runs=="1")
	ContR1<-subset(aggregate(R1[,4:6], by=list(R1[,3]), FUN=mean, na.rm=TRUE), Group.1=="cont")
	AveR1<-cbind(R1[,1:3], "group"=populations[I], "stage"="early", "0"= (((R1[,4]-ContR1[,2])/ContR1[,2])*100), "1"=(((R1[,5]-ContR1[,3])/ContR1[,3])*100), "5"=(((R1[,6]-ContR1[,4])/ContR1[,4])*100))
	R2<-subset(Ave1, Runs=="2" | Runs=="4")
	#R2<-subset(Ave1, Runs=="2")
	ContR2<-subset(aggregate(R2[,4:6], by=list(R2[,3]), FUN=mean, na.rm=TRUE), Group.1=="cont")
	AveR2<-cbind(R2[,1:3], "group"=populations[I], "stage"="late", "0"= (((R2[,4]-ContR2[,2])/ContR2[,2])*100), "1"=(((R2[,5]-ContR2[,3])/ContR2[,3])*100), "5"=(((R2[,6]-ContR2[,4])/ContR2[,4])*100))
	R3<-subset(Ave1, Runs=="3")
	ContR3<-subset(aggregate(R3[,4:6], by=list(R3[,3]), FUN=mean, na.rm=TRUE), Group.1=="cont")
	AveR3<-cbind(R3[,1:3], "group"=populations[I], "stage"="early", "0"= (((R3[,4]-ContR3[,2])/ContR3[,2])*100), "1"=(((R3[,5]-ContR3[,3])/ContR3[,3])*100), "5"=(((R3[,6]-ContR3[,4])/ContR3[,4])*100))
	R4<-subset(Ave1, Runs=="4")
	ContR4<-subset(aggregate(R4[,4:6], by=list(R4[,3]), FUN=mean, na.rm=TRUE), Group.1=="cont")
	AveR4<-cbind(R4[,1:3], "group"=populations[I], "stage"="late", "0"= (((R4[,4]-ContR4[,2])/ContR4[,2])*100), "1"=(((R4[,5]-ContR4[,3])/ContR4[,3])*100), "5"=(((R4[,6]-ContR4[,4])/ContR4[,4])*100))


	means[[I]]<- rbind(AveR1, AveR2)
}

Datta<-data.frame(do.call(rbind, means))
summary(Datta)


species<-c("madracis", "porites" , "siderastrea")
populations<-c("picophyto", "prochl" , "synPop1", "synPop2", "synPop3", "hetBact", "fullSyn")
collect<-list()
collected<-list()
for(G in 1:3)
{	spec<-subset(Datta, category==species[G])
	for(I in 1:7)
	{	group<-subset(spec, group==populations[I])
		tmSer<-wilcox.test((group[,8]+200) ~ group$stage, p.adjust="none", paired=FALSE)[[3]]
		ttest<-t.test((group[,8]+200) ~ group$stage)[[3]]
		shap<-shapiro.test((group[,8]+200))[[2]]
		active5<-na.omit(melt(group[,-2], id.vars=c("Sample","category", "group", "stage"), variable.name="time", value.name="change"))
		active5$Sample<-droplevels(active5$Sample)
		active5$category<-droplevels(active5$category)
		active5$group<-droplevels(active5$group)
		fit1 <- lmer(change ~ time*stage + (1|Sample), data=active5, REML=FALSE)
		tmSerT<-anova(fit1)[[6]][2]
		collect[[I]]<-cbind("spec"=species[G], "group"=populations[I], "wilcox"=tmSer, "ttest"=ttest, "shapiro"=shap, "mlm"=tmSerT)
	}
	collected[[G]]<-data.frame(do.call(rbind, collect))
}	
collected
write.csv(collected, "stats-figure5.csv")



me<-aggregate(Datta[,8], by=list(Datta$category, Datta$group, Datta$stage), FUN=mean, na.rm=TRUE)
se<-aggregate(Datta[,8], by=list(Datta$category, Datta$group, Datta$stage), FUN=sd, na.rm=TRUE)

madM<-subset(me, Group.1=="madracis")
madse<-subset(se, Group.1=="madracis")

porM<-subset(me, Group.1=="porites")
porse<-subset(se, Group.1=="porites")

sidM<-subset(me, Group.1=="siderastrea")
sidse<-subset(se, Group.1=="siderastrea")


c("#d9d9d9","#a1d99b","#969696","#41ab5d","#d9d9d9","#fd8d3c","#969696","#fc4e2a","#d9d9d9","#fd8d3c","#969696","#fc4e2a","#d9d9d9","#fd8d3c","#969696","#fc4e2a","#d9d9d9","#fd8d3c","#969696","#fc4e2a","#d9d9d9","#6baed6","#969696","#2171b5","#d9d9d9","#8c6bb1","#969696","#88419d")

COLOR<-list(c("#a1d99b","#41ab5d"),c("#6baed6","#2171b5"),c("#fd8d3c","#fc4e2a"),c("#fd8d3c","#fc4e2a"),c("#fd8d3c","#fc4e2a"),c("#fd8d3c","#fc4e2a"))

pdf(file = "early-vs-late.pdf", width =2.75, height = 2, bg="transparent")

par(mfrow=c(1,8), oma = c(0.1, 3, 0.1, 0.1), mar = c(0.75, 0.1, 0.75, 0.35))

populations<-c("picophyto", "prochl", "fullSyn","synPop1", "synPop2", "synPop3")
for(I in 1:6)
{
	means<-subset(madM, Group.2==populations[I])
	ses<-subset(madse, Group.2==populations[I])
	coll<-c(COLOR[[I]][1], COLOR[[I]][2])
	if(I==1)
	{
		bar<-barplot(means[,4], width=c(0.78),axes=FALSE, ylim=c(-100, 0), col=coll[means$Group.3])
		arrows(bar, ((means[,4])+ses[,4]/sqrt(1)), bar, ((means[,4])-ses[,4]/sqrt(1)), lwd=1.0, angle=90, code=3, length=0.025)
		axis(2, col.axis="black", las=2, cex.axis=1.1, at=c(seq(from=-90,to=0,by=45)))
		abline(h=0, lwd=1, col="black")
	}	
	if(I>1)
	{
		bar<-barplot(means[,4], width=c(0.78),axes=FALSE, ylim=c(-100, 0), col=coll[means$Group.3])
		arrows(bar, ((means[,4])+ses[,4]/sqrt(1)), bar, ((means[,4])-ses[,4]/sqrt(1)), lwd=1.0, angle=90, code=3, length=0.025)
		#axis(2, col.axis="black", las=2, cex.axis=1.1, at=c(seq(from=-90,to=0,by=30)))
	}

	box(col="black")
	grid(col="grey")
}

means<-subset(porM, Group.2==populations[1])
ses<-subset(porse, Group.2==populations[1])
coll<-c(COLOR[[1]][1], COLOR[[1]][2])
bar<-barplot(means[,4], width=c(0.78),axes=FALSE, ylim=c(-100, 0), col=coll[means$Group.3])
arrows(bar, ((means[,4])+ses[,4]/sqrt(1)), bar, ((means[,4])-ses[,4]/sqrt(1)), lwd=1.0, angle=90, code=3, length=0.025)
#axis(2, col.axis="black", las=2, cex.axis=1.1, at=c(seq(from=-90,to=0,by=30)))
box(col="black")
grid(col="grey")

means<-subset(porM, Group.2==populations[3])
ses<-subset(porse, Group.2==populations[3])
coll<-c(COLOR[[3]][1], COLOR[[3]][2])
bar<-barplot(means[,4], width=c(0.78),axes=FALSE, ylim=c(-100, 0), col=coll[means$Group.3])
arrows(bar, ((means[,4])+ses[,4]/sqrt(1)), bar, ((means[,4])-ses[,4]/sqrt(1)), lwd=1.0, angle=90, code=3, length=0.025)
#axis(2, col.axis="black", las=2, cex.axis=1.1, at=c(seq(from=-90,to=0,by=30)))
box(col="black")
grid(col="grey")

dev.off()









