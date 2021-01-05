library(geiger); library(geomorph)
setwd("/Users/Home/Dropbox/Amherst/Project/Moduarity/MarsupialPlacental/3D/Data/Carnivores/ReAnalyze/")

AllCarns<-readland.tps("Individual-MarPlaJaw_20200805.tps", specID = "ID")

####Procrustes####
Y.gpa<-gpagen(A = AllCarns)

####Allometry####
gdf <- geomorph.data.frame(Y.gpa, species = dimnames(Y.gpa$coords)[[3]])

##Obtaining size-adjusted residuals (and allometry-free shapes)##
HybAnova <- procD.lm(coords~Csize, data = gdf, iter = 999, RRPP=TRUE) 
summary(HybAnova) 

shape.resid <- arrayspecs(HybAnova$residuals, p=dim(Y.gpa$coords)[1], k=dim(Y.gpa$coords)[2]) # size-adjusted residuals

#allometry-free shapes
MandibleShape <- shape.resid + array(Y.gpa$consensus, dim(shape.resid))


##Save PC Scores
#Read in Grouping data and TPS file
CarnData<-read.csv("Individual-Data_20200829.csv")

#Match TPS order to data file order
MandibleShape<-MandibleShape[,,(match(CarnData$SpecimenID, dimnames(MandibleShape)[[3]]))]
MandiblePCA<-gm.prcomp(A = MandibleShape)

#Assign colors
CarnData$GroupB<-as.factor(CarnData$GroupB)
COLVEC<-c("#8da0cb","#fc8d62","#66c2a5")

par(pty='s')
plot(MandiblePCA, axis1 = 1, axis2 = 2, pch=1, col=COLVEC[unclass(CarnData$GroupB)])
text(MandiblePCA$x, pos = 4, label = dimnames(MandibleShape)[[3]])

write.csv(MandiblePCA$x, "ManPCScores.csv")


####Mean TPS file####
#Coordinates#
DataNak <-read.csv("Individual-Data_20200829.csv")
DataNak$PhyloName<-as.factor(DataNak$PhyloName)

#Create mean GPA from all specimens
CarnMeanShapes<-vector("list", length(levels(DataNak$PhyloName)))

for(i in 1:length(levels(DataNak$PhyloName))){
  
  CarnNames <- levels(DataNak$PhyloName)[i]
  NamesRows <- which(DataNak$PhyloName == CarnNames)
  
  CatTPS <- MandibleShape[,,as.character(DataNak[NamesRows,2])]
  
  #If you have a single entry need to convert to array (R thinks it is a matrix)
  if(class(CatTPS) != "array"){
    CatTPS<-array(CatTPS, dim=c(12,3,1))
  }
  
  Carn.gpa<-gpagen(CatTPS)
  CarnMeanShapes[[i]] <- mshape(Carn.gpa$coords)
  
}

#Convert from list to array
CarnMeanShapes <-array(unlist(CarnMeanShapes), dim=c(12,3,length(levels(DataNak$PhyloName))))

#Add species names
CarnDimName<-vector(mode="list", 3)
CarnDimName[[3]]<-levels(DataNak$PhyloName)
dimnames(CarnMeanShapes)<-CarnDimName

writeland.tps(CarnMeanShapes, "Mean-MarPla-Shape_20200816.tps")


##Save PC Scores
#Read in Grouping data and TPS file
CarnMeanShapes<-readland.tps(file = "Mean-MarPla-Shape_20200816.tps", specID = "ID")
CarnMeanData<-read.csv("Mean-Data_20200806.csv")

#Match TPS order to data file order
CarnMeanShapes<-CarnMeanShapes[,,(match(CarnMeanData$PhyloID, dimnames(CarnMeanShapes)[[3]]))]
MandibleMeanPCA<-gm.prcomp(A = CarnMeanShapes)

#Assign colors
CarnMeanData$GroupB<-as.factor(CarnMeanData$GroupB)
COLVEC<-c("#8da0cb","#fc8d62","#66c2a5")

par(pty='s')
plot(MandibleMeanPCA, axis1 = 1, axis2 = 2, pch=1, col=COLVEC[unclass(CarnMeanData$GroupB)])
text(MandibleMeanPCA$x, pos = 4, label = dimnames(CarnMeanShapes)[[3]])

write.csv(MandibleMeanPCA$x, "ManMeanPCScores.csv")



####Tree data engineering####
setwd("/Users/Home/Dropbox/Amherst/Project/Moduarity/MarsupialPlacental/3D/Data/Carnivores/ReAnalyze/")
#Read distribution of trees
UltraMammal<-read.nexus(file = "UltraMammal1K_2019.nex")
CarnMeanData<-read.csv("Mean-Data_20200806.csv")
CarnMeanShapes<-readland.tps(file = "Mean-MarPla-Shape_20200816.tps", specID = "ID")

#Retain taxa in tree
CarnMeanShapesTr<-CarnMeanShapes[,,dimnames(CarnMeanShapes)[[3]][dimnames(CarnMeanShapes)[[3]]%in%UltraMammal$tree_4444$tip.label]]

#Match TPS order to tree order
CarnMeanShapesTr<-CarnMeanShapesTr[,,(match(UltraMammal$tree_4444$tip.label, dimnames(CarnMeanShapesTr)[[3]]))]



####Disparity####
#Setup for disparity analysis
gp.end <- as.factor(CarnMeanData$GroupA)
gdf <- geomorph.data.frame(coords=CarnMeanShapesTr, gp.end = gp.end)

Mammalpgls = lapply(UltraMammal, function(x) {
  procD.pgls(f1 = coords~1, phy = x, data = gdf, iter = 999) # phylogenetic generalized least squares
})

Mammaldisparity.GpA = lapply(Mammalpgls, function(x) {
  morphol.disparity(f1 = x, groups = ~ gp.end, data = gdf, iter = 999)
})

#Collate Data
save(Mammalpgls,file="Mammalpgls.rda")
save(Mammaldisparity.GpA,file="Mammaldisparity_GpA.rda")
save(Mammaldisparity.GpB,file="Mammaldisparity_GpB.rda")

#PValue
PVDistPval.MatrixMean<-list()
for(i in 1:1000){
  PVDistPval.MatrixMean[[i]]<-Mammaldisparity.GpA[[i]]$PV.dist.Pval
}
apply(simplify2array(PVDistPval.MatrixMean), 1:2, mean)
apply(simplify2array(PVDistPval.MatrixMean), 1:2, median)
apply(simplify2array(PVDistPval.MatrixMean), 1:2, sd)
apply(simplify2array(PVDistPval.MatrixMean), 1:2, quantile, probs=c(0.025,0.975))

#Procrustes Variances
PVar.MatrixMean<-list()
for(i in 1:1000){
  PVar.MatrixMean[[i]]<-Mammaldisparity.GpA[[i]]$Procrustes.var
}
PVar.Mean<-matrix(data = unlist(PVar.MatrixMean), nrow = 1000, ncol = 2, byrow = T)

apply(PVar.Mean, 2, mean)
apply(PVar.Mean, 2, quantile, probs=c(0.025,0.975))
write.csv(PVar.Mean, "PVarMean_CarDas-1K.csv")

###Plot Disparity
library(ggplot2)
Mamdi<-read.csv("Results/Disparity/PVarMean_CanDasFel-1K.csv")
Mamdi<-read.csv("Results/Disparity/PVarMean_CarDas-1K.csv")

##Car-Das
#Violin plots
ggplot(Mamdi, aes(x=factor(Group), y=PV, fill=factor(Group))) + geom_violin(width=1, trim=F) + labs(title="Disparity", y = "Procrustes Variances", x = "Clade", fill = "Clade") + theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "gray"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "gray"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_fill_discrete(name="Clade", breaks=c("Das", "Ecar"), labels=c("Dasyuromorphia", "Carnivora"))

#Boxplot
ggplot(Mamdi, aes(x=factor(Group), y=PV, fill=factor(Group))) + geom_boxplot() + labs(title="Disparity", y = "Procrustes Variances", x = "Clade", fill = "Clade") + theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "gray"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "gray"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_fill_discrete(name="Clade", breaks=c("Das", "Ecar"), labels=c("Dasyuromorphia", "Carnivora"))

##Das-Can-Fel
#Violin plots
ggplot(Mamdi, aes(x=factor(Group), y=PV, fill=factor(Group))) + geom_violin(width=1, trim=F) + labs(title="Disparity", y = "Procrustes Variances", x = "Clade", fill = "Clade") + theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "gray"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "gray"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_fill_discrete(name="Clade", breaks=c("Das", "Ecan", "Fel"), labels=c("Dasyuromorphia", "Caniformia", "Feliformia"))

#Boxplot
ggplot(Mamdi, aes(x=factor(Group), y=PV, fill=factor(Group))) + geom_boxplot() + labs(title="Disparity", y = "Procrustes Variances", x = "Clade", fill = "Clade") + theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "gray"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "gray"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_fill_discrete(name="Clade", breaks=c("Das", "Ecan", "Fel"), labels=c("Dasyuromorphia", "Caniformia", "Feliformia"))


####Rates####
#Rates between clades
gp.end <- as.factor(CarnMeanData$GroupB)
names(gp.end) <- CarnMeanData$PhyloID

Mammalcompareratesmulti.GpB = lapply(UltraMammal, function(x) {
  compare.evol.rates(CarnMeanShapesTr, x, gp.end, iter = 999, method = "simulation", print.progress = TRUE)
})

#save(Mammalcompareratesmulti.GpB,file="Mammalrates_GpB.rda")

##Collate Rate Data Group B
#Sigma D
Rate.GpB<-matrix(data = NA, nrow = 1000, ncol = 3)
for(i in 1:1000){
  Rate.GpB[i,]<-as.numeric(Mammalcompareratesmulti.GpB[[i]]$sigma.d.gp)
}
colMeans(Rate.GpB)
apply(X = Rate.GpB, MARGIN = 2, FUN = quantile, probs=c(0.025,0.975))
write.csv(Rate.GpB, "RatesOutput_GpB.csv")

#P-Value
Rate.GpB.P<-matrix(data = NA, nrow = 1000, ncol = 3)
for(i in 1:1000){
  Rate.GpB.P[i,]<-as.numeric(Mammalcompareratesmulti.GpB[[i]]$pairwise.pvalue[c(1,2,3)])
}
colMeans(Rate.GpB.P)
apply(X = Rate.GpB.P, MARGIN = 2, FUN = quantile, probs=c(0.025,0.975))
write.csv(Rate.GpB.P, "RatesOutputP_GpB.csv")

#GpA
gp.end <- as.factor(CarnMeanData$GroupA)
names(gp.end) <- CarnMeanData$PhyloID

Mammalcompareratesmulti.GpA = lapply(UltraMammal, function(x) {
  compare.evol.rates(CarnMeanShapesTr, x, gp.end, iter = 999, method = "simulation", print.progress = TRUE)
})

#save(Mammalcompareratesmulti.GpA,file="Mammalrates_GpA.rda")

##Collate Rate Data Group A
#Sigma D
Rate.GpA<-matrix(data = NA, nrow = 1000, ncol = 2)
for(i in 1:1000){
  Rate.GpA[i,]<-as.numeric(Mammalcompareratesmulti.GpA[[i]]$sigma.d.gp)
}
colMeans(Rate.GpA)
apply(X = Rate.GpA, MARGIN = 2, FUN = quantile, probs=c(0.025,0.975))
write.csv(Rate.GpA, "RatesOutput_GpA.csv")

#P-Value
Rate.GpA.P<-NULL
for(i in 1:1000){
  Rate.GpA.P[i]<-as.numeric(Mammalcompareratesmulti.GpA[[i]]$pairwise.pvalue)
}
mean(Rate.GpA.P)
quantile(Rate.GpA.P, probs=c(0.025,0.975))
write.csv(Rate.GpA.P, "RatesOutputP_GpA.csv")

###Plot Rates
library(ggplot2)
Mamra<-read.csv("Results/Rates/RatesMean_CanDasFel-1K.csv")
Mamra<-read.csv("ResultsNew/Rates/RatesMean_CarDas-1K_New.csv")

##Car-Das
#Violin plots
ggplot(Mamra, aes(x=factor(Group), y=Rate, fill=factor(Group))) + geom_violin(width=1, trim=F) + labs(title="Rates", y = bquote('Brownian rate ('*~sigma^2*')'), x = "Clade", fill = "Clade") + theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "gray"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "gray"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_fill_discrete(name="Clade", breaks=c("Das", "Ecan"), labels=c("Dasyuromorphia", "Carnivora"))

#Boxplot
ggplot(Mamra, aes(x=factor(Group), y=Rate, fill=factor(Group))) + geom_boxplot() + labs(title="Rates", y = bquote('Brownian rate ('*~sigma^2*')'), x = "Clade", fill = "Clade") + theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "gray"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "gray"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_fill_discrete(name="Clade", breaks=c("Das", "Ecan"), labels=c("Dasyuromorphia", "Carnivora"))


##Das-Can-Fel
#Violin plots
ggplot(Mamra, aes(x=factor(Group), y=Rate, fill=factor(Group))) + geom_violin(width=1, trim=F) + labs(title="Rates", y = bquote('Brownian rate ('*~sigma^2*')'), x = "Clade", fill = "Clade") + theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "gray"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "gray"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_fill_discrete(name="Clade", breaks=c("Das", "Ecan", "Fel"), labels=c("Dasyuromorphia", "Caniformia", "Feliformia"))

#Boxplot
ggplot(Mamra, aes(x=factor(Group), y=Rate, fill=factor(Group))) + geom_boxplot() + labs(title="Rates", y = bquote('Brownian rate ('*~sigma^2*')'), x = "Clade", fill = "Clade") + theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "gray"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "gray"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_fill_discrete(name="Clade", breaks=c("Das", "Ecan", "Fel"), labels=c("Dasyuromorphia", "Caniformia", "Feliformia"))


####EMMLi Phylo####
library(EMMLi); library(paleotree); library(paleomorph)
UltraMammal<-read.nexus(file = "UltraMammal1K_2019.nex")
CarnMeanData<-read.csv("Mean-Data_20200806.csv")
CarnMeanShapes<-readland.tps(file = "Mean-MarPlaShape_20200805.tps", specID = "ID")

#Retain taxa in tree
CarnMeanShapesTr<-CarnMeanShapes[,,dimnames(CarnMeanShapes)[[3]][dimnames(CarnMeanShapes)[[3]]%in%UltraMammal$tree_6920$tip.label]]

#Match TPS order to tree order
CarnMeanShapesTr<-CarnMeanShapesTr[,,(match(UltraMammal$tree_6920$tip.label, dimnames(CarnMeanShapesTr)[[3]]))]

#restart here
CarnMeanShapesTrEMMLi<-CarnMeanShapesTr

#Subset by clades
KeepDas<-CarnMeanData[CarnMeanData$GroupB=="Das",]
KeepCar<-CarnMeanData[CarnMeanData$GroupA=="Car",]
KeepCan<-CarnMeanData[CarnMeanData$GroupB=="Can",]
KeepFel<-CarnMeanData[CarnMeanData$GroupB=="Fel",]

CarnMeanShapesTrEMMLi<-CarnMeanShapesTrEMMLi[,,KeepDas$PhyloID]
CarnMeanShapesTrEMMLi<-CarnMeanShapesTrEMMLi[,,KeepCar$PhyloID]
CarnMeanShapesTrEMMLi<-CarnMeanShapesTrEMMLi[,,KeepCan$PhyloID]
CarnMeanShapesTrEMMLi<-CarnMeanShapesTrEMMLi[,,KeepFel$PhyloID]

#Prune Trees
EMMLiDropTree<-lapply(UltraMammal,keep.tip,tip=dimnames(CarnMeanShapesTrEMMLi)[[3]])
class(EMMLiDropTree)<-"multiPhylo"


#Start the EMMLi analysis
coords.2d<-two.d.array(CarnMeanShapesTrEMMLi)

coords.contrasts = lapply(EMMLiDropTree, function(x) {
  apply(coords.2d, 2, FUN=pic, phy=x)
})

contrasts3d = lapply(coords.contrasts, function(x) {
  arrayspecs(x, p=dim(CarnMeanShapesTrEMMLi)[1], k = 3)
})

cormat.contrasts = lapply(contrasts3d, function(x) {
  dotcorr(x)
})


LMList<-read.csv("EMMLi/LMsmodel.csv")
LMList[,1]<-as.character(LMList[,1])

##With n-1
#Das, 33; Car, 206; Can, 108; Fel, 98.
Phylo.emmli.Das = lapply(cormat.contrasts, function(x) {
  EMMLi(corr = x, N_sample = 32, mod = LMList, pprob = 0.00000000001)
})

Phylo.emmli.Car = lapply(cormat.contrasts, function(x) {
  EMMLi(corr = x, N_sample = 205, mod = LMList, pprob = 0.00000000001)
})

Phylo.emmli.Can = lapply(cormat.contrasts, function(x) {
  EMMLi(corr = x, N_sample = 107, mod = LMList, pprob = 0.00000000001)
})

Phylo.emmli.Fel = lapply(cormat.contrasts, function(x) {
  EMMLi(corr = x, N_sample = 97, mod = LMList, pprob = 0.00000000001)
})


#Save EMMLi
save(Phylo.emmli.Das, file = "UltrametricTroph_PhyloemmliDas.rda")
save(Phylo.emmli.Car, file = "UltrametricTroph_PhyloemmliCar.rda")
save(Phylo.emmli.Can, file = "UltrametricTroph_PhyloemmliCan.rda")
save(Phylo.emmli.Fel, file = "UltrametricTroph_PhyloemmliFel.rda")


####Mean EMMLi Table####
##Produce mean summary table
library(plyr)
#Dasyuromorphia
ojEMMLiDas<-list()
for(i in 1:1000){
  ojEMMLiDas[[i]]<-Phylo.emmli.Das[[i]]$results
  rownames(ojEMMLiDas[[i]]) <- c()
}

EffectMeanDas <- as.data.frame(aaply(laply(ojEMMLiDas, as.matrix), c(2, 3), mean))
#EffectQuantDas <- aaply(laply(ojEMMLiDas, as.matrix), c(2, 3), quantile, probs=c(0.025,0.975))
rownames(EffectMeanDas)<-rownames(Phylo.emmli.Das[[1]]$results)
write.csv(EffectMeanDas,"EMMLiDasResultsTableMean.csv")

#Carnivora
ojEMMLiCar<-list()
for(i in 1:1000){
  ojEMMLiCar[[i]]<-Phylo.emmli.Car[[i]]$results
  rownames(ojEMMLiCar[[i]]) <- c()
}

EffectMeanCar <- as.data.frame(aaply(laply(ojEMMLiCar, as.matrix), c(2, 3), mean))
#EffectQuantCar <- aaply(laply(ojEMMLiCar, as.matrix), c(2, 3), quantile, probs=c(0.025,0.975))
rownames(EffectMeanCar)<-rownames(Phylo.emmli.Car[[1]]$results)
write.csv(EffectMeanCar,"EMMLiCarResultsTableMean.csv")

#Caniformia
ojEMMLiCan<-list()
for(i in 1:1000){
  ojEMMLiCan[[i]]<-Phylo.emmli.Can[[i]]$results
  rownames(ojEMMLiCan[[i]]) <- c()
}

EffectMeanCan <- as.data.frame(aaply(laply(ojEMMLiCan, as.matrix), c(2, 3), mean))
#EffectQuantCan <- aaply(laply(ojEMMLiCan, as.matrix), c(2, 3), quantile, probs=c(0.025,0.975))
rownames(EffectMeanCan)<-rownames(Phylo.emmli.Can[[1]]$results)
write.csv(EffectMeanCan,"EMMLiCanResultsTableMean.csv")

#Feliformia
ojEMMLiFel<-list()
for(i in 1:1000){
  ojEMMLiFel[[i]]<-Phylo.emmli.Fel[[i]]$results
  rownames(ojEMMLiFel[[i]]) <- c()
}

EffectMeanFel <- as.data.frame(aaply(laply(ojEMMLiFel, as.matrix), c(2, 3), mean))
#EffectQuantFel <- aaply(laply(ojEMMLiFel, as.matrix), c(2, 3), quantile, probs=c(0.025,0.975))
rownames(EffectMeanFel)<-rownames(Phylo.emmli.Fel[[1]]$results)
write.csv(EffectMeanFel,"EMMLiFelResultsTableMean.csv")


#######EMMLi module consistency####
##Das
##Zero AICc units
EMsmodulesZero<-list()
for(i in 1:1000){
  EMsmodulesZero[[i]]<-names(which(Phylo.emmli.Das[[i]]$results[,5]==0))
}
EMsmodulesZeroEdit = lapply(EMsmodulesZero, function(x) {
  unique(gsub("^([^.]+.[^.]+).*","\\1",x))
})
table(unlist(EMsmodulesZeroEdit))

##Two AICc units
EMsmodulesTwo<-list()
for(i in 1:1000){
  EMsmodulesTwo[[i]]<-names(which(Phylo.emmli.Das[[i]]$results[,5]<2))
}
EMsmodulesTwoEdit = lapply(EMsmodulesTwo, function(x) {
  unique(gsub("^([^.]+.[^.]+).*","\\1",x))
})
table(unlist(EMsmodulesTwoEdit))

##Car
##Zero AICc units
EMsmodulesZero<-list()
for(i in 1:1000){
  EMsmodulesZero[[i]]<-names(which(Phylo.emmli.Car[[i]]$results[,5]==0))
}
EMsmodulesZeroEdit = lapply(EMsmodulesZero, function(x) {
  unique(gsub("^([^.]+.[^.]+).*","\\1",x))
})
table(unlist(EMsmodulesZeroEdit))

##Two AICc units
EMsmodulesTwo<-list()
for(i in 1:1000){
  EMsmodulesTwo[[i]]<-names(which(Phylo.emmli.Car[[i]]$results[,5]<2))
}
EMsmodulesTwoEdit = lapply(EMsmodulesTwo, function(x) {
  unique(gsub("^([^.]+.[^.]+).*","\\1",x))
})
table(unlist(EMsmodulesTwoEdit))


##Can
##Zero AICc units
EMsmodulesZero<-list()
for(i in 1:1000){
  EMsmodulesZero[[i]]<-names(which(Phylo.emmli.Can[[i]]$results[,5]==0))
}
EMsmodulesZeroEdit = lapply(EMsmodulesZero, function(x) {
  unique(gsub("^([^.]+.[^.]+).*","\\1",x))
})
table(unlist(EMsmodulesZeroEdit))

##Two AICc units
EMsmodulesTwo<-list()
for(i in 1:1000){
  EMsmodulesTwo[[i]]<-names(which(Phylo.emmli.Can[[i]]$results[,5]<2))
}
EMsmodulesTwoEdit = lapply(EMsmodulesTwo, function(x) {
  unique(gsub("^([^.]+.[^.]+).*","\\1",x))
})
table(unlist(EMsmodulesTwoEdit))


##Fel
##Zero AICc units
EMsmodulesZero<-list()
for(i in 1:1000){
  EMsmodulesZero[[i]]<-names(which(Phylo.emmli.Fel[[i]]$results[,5]==0))
}
EMsmodulesZeroEdit = lapply(EMsmodulesZero, function(x) {
  unique(gsub("^([^.]+.[^.]+).*","\\1",x))
})
table(unlist(EMsmodulesZeroEdit))

##Two AICc units
EMsmodulesTwo<-list()
for(i in 1:1000){
  EMsmodulesTwo[[i]]<-names(which(Phylo.emmli.Fel[[i]]$results[,5]<2))
}
EMsmodulesTwoEdit = lapply(EMsmodulesTwo, function(x) {
  unique(gsub("^([^.]+.[^.]+).*","\\1",x))
})
table(unlist(EMsmodulesTwoEdit))


##Plotting
CanDasEm<-read.csv("Results/EMMLi/Mam_EMMLiSummaryPlot.csv")

# Stacked + percent
ggplot(CanDasEm, aes(fill=Module, y=Frequency, x=AICc)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw() + facet_wrap(~Clade)
