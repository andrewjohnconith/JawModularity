library(geiger); library(geomorph); library(paleomorph); library(EMMLi); library(phytools); library(plotrix)
setwd("/Users/home/Dropbox/Amherst/Project/Moduarity/MarsupialPlacental/3D/Paper/Submission/ResubResults/Writing/Final/GitHub/")

MandibleShape<-readland.tps("Individual-MarPlaJaw_20200805.tps", specID = "ID")
DataNak <-read.csv("Individual-Data_20200829_order.csv")
morphmean <- aggregate(CentroidSize~PhyloName,data=DataNak,FUN="mean",na.rm=TRUE,na.action=NULL)

Y.gpa<-gpagen(A = MandibleShape)
ManPCA<-gm.prcomp(Y.gpa$coords)

COLVEC<-c("#a6cee3", "#b2df8a", "#1f78b4")
DataNak$GroupB<-as.factor(DataNak$GroupB)

par(pty='s') #Makes plot square
pc.plot<-plot(ManPCA, pch=19, col=COLVEC[unclass(DataNak$GroupB)])
shapeHulls(pc.plot, groups = DataNak$GroupB, group.cols = c("#a6cee3", "#1f78b4", "#b2df8a"))
legend("topleft", legend = levels(DataNak$GroupB), col = COLVEC, pch = 19, horiz = T)
#text(ManPCA$x[,1],ManPCA$x[,2],labels = DataNak$GroupB, pos=4)

####Mean TPS file####
#Coordinates#
DataNak <-read.csv("Individual-Data_20200829_order.csv")
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


####Tree data engineering####
#Read distribution of trees
UltraMammal<-read.nexus(file = "UltraMammal1K_2019.nex")
CarnMeanData<-read.csv("Mean-Data_20210318.csv")
CarnMeanData$GroupA<-as.factor(CarnMeanData$GroupA)
CarnMeanData$GroupB<-as.factor(CarnMeanData$GroupB)

Y.gpa<-gpagen(A = CarnMeanShapes)
CarnMeanShapes<-Y.gpa$coords

#Retain taxa in tree
CarnMeanShapesTr<-CarnMeanShapes[,,dimnames(CarnMeanShapes)[[3]][dimnames(CarnMeanShapes)[[3]]%in%UltraMammal$tree_4444$tip.label]]

#Match TPS order to tree order
CarnMeanShapesTr<-CarnMeanShapesTr[,,(match(UltraMammal$tree_4444$tip.label, dimnames(CarnMeanShapesTr)[[3]]))]

#Procrustes and add back in mean centroid size
CarnMeanShapesTr<-gpagen(CarnMeanShapesTr)
CarnMeanShapesTr$Csize<-CarnMeanData$MeanCentroidSize

#Phylogenetic ANOVA
MyPhyANOVA_A<-phylANOVA(tree = UltraMammal$tree_4444, x = CarnMeanData$GroupA, y = CarnMeanData$MeanCentroidSize, nsim = 999)
boxplot(CarnMeanData$MeanCentroidSize~CarnMeanData$GroupA)
MyPhyANOVA_B<-phylANOVA(tree = UltraMammal$tree_4444, x = CarnMeanData$GroupB, y = CarnMeanData$MeanCentroidSize, posthoc = T, p.adj = "holm", nsim = 999)
boxplot(CarnMeanData$MeanCentroidSize~CarnMeanData$GroupB)

##plot mean shapes
#ManPCA<-gm.prcomp(A = CarnMeanShapesTr$coords, phy = UltraMammal$tree_4444, align.to.phy = T)
ManPCA<-gm.prcomp(A = CarnMeanShapesTr$coords, phy = UltraMammal$tree_4444, align.to.phy = F)

COLVEC<-c("#a6cee3", "#b2df8a", "#1f78b4")

par(pty='s') #Makes plot square
pc.plot<-plot(ManPCA, pch=19, col=COLVEC[unclass(CarnMeanData$GroupB)])
shapeHulls(pc.plot, groups = CarnMeanData$GroupB, group.cols = c("#a6cee3", "#1f78b4", "#b2df8a"))
legend("bottomleft", legend = levels(CarnMeanData$GroupB), col = COLVEC, pch = 19, horiz = T)
text(ManPCA$x[,1],ManPCA$x[,2],labels = CarnMeanData$PhyloID, pos=4)



####Allometry####
#Setup for allometry analysis
gp.end <- as.factor(CarnMeanData$GroupA)
gdf <- geomorph.data.frame(CarnMeanShapesTr, gp.end = gp.end)

fit.size = lapply(UltraMammal, function(x) {
  procD.pgls(coords~log(Csize), phy = x,
             data = gdf, print.progress = F)
})


####Disparity####
##Setup for disparity analysis
#Group A
gp.end <- as.factor(CarnMeanData$GroupA)
gdf <- geomorph.data.frame(CarnMeanShapesTr, gp.end = gp.end)


Mammalpgls.GpA = lapply(UltraMammal, function(x) {
  procD.pgls(f1 = coords~log(Csize)+gp.end, phy = x, data = gdf, iter = 999, print.progress = T) # phylogenetic generalized least squares
})


Mammaldisparity.GpA = lapply(Mammalpgls.GpA, function(x) {
  morphol.disparity(f1 = x, groups = ~ gp.end, data = gdf, iter = 999, print.progress = T)
})

#Collate Disparity Data
##P-value
Disparity_P<-as.data.frame(matrix(unlist(lapply(Mammaldisparity.GpA, "[[", 2)), nrow = 1000, ncol = 4, byrow = T))[,2]

#Mean loop
PVDistPval.MatrixMean<-list()
for(i in 1:1000){
  PVDistPval.MatrixMean[[i]]<-Mammaldisparity.GpA[[i]]$PV.dist.Pval
}

apply(simplify2array(PVDistPval.MatrixMean), 1:2, mean)
apply(simplify2array(PVDistPval.MatrixMean), 1:2, median)
apply(simplify2array(PVDistPval.MatrixMean), 1:2, sd)
apply(simplify2array(PVDistPval.MatrixMean), 1:2, quantile, probs=c(0.025,0.975))

##Procrustes Variances
Disparity_ProcVar<-as.data.frame(matrix(unlist(lapply(Mammaldisparity.GpA, "[[", 1)), nrow = 1000, ncol = 2, byrow = T))

#Mean Loop
PVar.MatrixMean<-list()
for(i in 1:1000){
  PVar.MatrixMean[[i]]<-Mammaldisparity.GpA[[i]]$Procrustes.var
}
PVar.Mean<-matrix(data = unlist(PVar.MatrixMean), nrow = 1000, ncol = 2, byrow = T)

apply(PVar.Mean, 2, mean)
apply(PVar.Mean, 2, quantile, probs=c(0.025,0.975))

#Collate
Disparity_Data<-cbind.data.frame(Disparity_ProcVar, Disparity_P)
colnames(Disparity_Data)<-c("Car_ProcV","Das_ProcV","CanDas_P")
write.csv(Disparity_Data, "Disparity_Data_CarDas-1K.csv")


#Group B
gp.end <- as.factor(CarnMeanData$GroupB)
gdf <- geomorph.data.frame(CarnMeanShapesTr, gp.end = gp.end)


Mammalpgls.GpB = lapply(UltraMammal, function(x) {
  procD.pgls(f1 = coords~Csize+gp.end, phy = x, data = gdf, iter = 999, print.progress = T) # phylogenetic generalized least squares
})

#Read in fit common group B
Mammaldisparity.GpB = lapply(Mammalpgls.GpB, function(x) {
  morphol.disparity(f1 = x, groups = ~ gp.end, data = gdf, iter = 999, print.progress = T)
})


##P-value
#Can-Das, Can-Fel, Das-Fel
Disparity_P<-as.data.frame(matrix(unlist(lapply(Mammaldisparity.GpB, "[[", 2)), nrow = 1000, ncol = 9, byrow = T))[,c(2,3,6)]

#Mean loop (Does the same thing as the above line)
PVDistPval.MatrixMean<-list()
for(i in 1:1000){
  PVDistPval.MatrixMean[[i]]<-Mammaldisparity.GpB[[i]]$PV.dist.Pval
}

apply(simplify2array(PVDistPval.MatrixMean), 1:2, mean)
apply(simplify2array(PVDistPval.MatrixMean), 1:2, median)
apply(simplify2array(PVDistPval.MatrixMean), 1:2, sd)
apply(simplify2array(PVDistPval.MatrixMean), 1:2, quantile, probs=c(0.025,0.975))

##Procrustes Variances
Disparity_ProcVar<-as.data.frame(matrix(unlist(lapply(Mammaldisparity.GpB, "[[", 1)), nrow = 1000, ncol = 3, byrow = T))

#Mean Loop
PVar.MatrixMean<-list()
for(i in 1:1000){
  PVar.MatrixMean[[i]]<-Mammaldisparity.GpB[[i]]$Procrustes.var
}
PVar.Mean<-matrix(data = unlist(PVar.MatrixMean), nrow = 1000, ncol = 3, byrow = T)

apply(PVar.Mean, 2, mean)
apply(PVar.Mean, 2, quantile, probs=c(0.025,0.975))

#Collate
Disparity_Data<-cbind.data.frame(Disparity_ProcVar, Disparity_P)
colnames(Disparity_Data)<-c("Can_ProcV","Das_ProcV","Fel_ProcV","CanDas_P", "CanFel_P", "DasFel_P")
write.csv(Disparity_Data, "Disparity_Data_CanFelDas-1K.csv")


##Compare disparity metrics
setwd("/Users/home/Dropbox/Amherst/Project/Moduarity/MarsupialPlacental/3D/Paper/Submission/ResubResults/Writing/Final/GitHub/")
#setwd("/Users/Home/Dropbox/Amherst/Project/Moduarity/MarsupialPlacental/3D/Paper/Submission/ResubResults/Disparity/")
DisparityA<-read.csv("Disparity_Data_CarDas-1K.csv", row.names = 1)
colMeans(DisparityA); apply(X = DisparityA, MARGIN = 2, FUN = median, na.rm = T)

apply(X = DisparityA, MARGIN = 2, FUN = std.error, na.rm = T)
apply(X = DisparityA, MARGIN = 2, FUN = quantile, probs=c(0.025,0.975), na.rm = T)

DisparityB<-read.csv("Disparity_Data_CanFelDas-1K.csv", row.names = 1)
colMeans(DisparityB); apply(X = DisparityB, MARGIN = 2, FUN = median, na.rm = T)

apply(X = DisparityB, MARGIN = 2, FUN = std.error, na.rm = T)
apply(X = DisparityB, MARGIN = 2, FUN = quantile, probs=c(0.025,0.975), na.rm = T)


###Plot Disparity
library(ggplot2)
setwd("/Users/home/Dropbox/Amherst/Project/Moduarity/MarsupialPlacental/3D/Paper/Submission/ResubResults")

##Car-Das
Mamdi<-read.csv("Disparity/Disparity_Data_CarDas-1Kstack.csv")
#Violin plots
ggplot(Mamdi, aes(x=factor(Group), y=PV, fill=factor(Group))) + geom_violin(width=1, trim=F) + labs(title="Disparity", y = "Procrustes Variances", x = "Clade", fill = "Clade") + theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "gray"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "gray"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_fill_discrete(name="Clade", breaks=c("Das", "Ecar"), labels=c("Dasyuromorphia", "Carnivora"))

##Das-Can-Fel
Mamdi<-read.csv("Disparity/Disparity_Data_CanFelDas-1Kstack.csv")
#Violin plots
ggplot(Mamdi, aes(x=factor(Group), y=PV, fill=factor(Group))) + geom_violin(width=1, trim=F) + labs(title="Disparity", y = "Procrustes Variances", x = "Clade", fill = "Clade") + theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "gray"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "gray"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_fill_discrete(name="Clade", breaks=c("Das", "Ecan", "Fel"), labels=c("Dasyuromorphia", "Caniformia", "Feliformia"))


####Rates####
#Rates between clades - GpA
gp.end <- as.factor(CarnMeanData$GroupA)
names(gp.end) <- CarnMeanData$PhyloID

Mammalcompareratesmulti.GpA = lapply(UltraMammal, function(x) {
  compare.evol.rates(A = CarnMeanShapesTr$coords, phy = x, gp = gp.end, iter = 999, method = "simulation", print.progress = T)
})

Rate_Ratio<-unlist(lapply(Mammalcompareratesmulti.GpA, "[[", "sigma.d.ratio"))
Rate_P<-unlist(lapply(Mammalcompareratesmulti.GpA, "[[", "P.value"))
Rate_ES<-unlist(lapply(Mammalcompareratesmulti.GpA, "[[", "Z"))
Rate_Clade<-as.data.frame(matrix(unlist(lapply(Mammalcompareratesmulti.GpA, "[[", "sigma.d.gp")), nrow = 1000, ncol = 2, byrow = T))

#Collate
Rate_Data<-cbind.data.frame(Rate_Ratio,Rate_P,Rate_ES,Rate_Clade)
colnames(Rate_Data)<-c("sigmadratio","P","Z","CarSigma", "DasSigma")
write.csv(Rate_Data, "Rate_Data_CarDas-1K.csv")


#Rates between clades - GpB
gp.end <- as.factor(CarnMeanData$GroupB)
names(gp.end) <- CarnMeanData$PhyloID

Mammalcompareratesmulti.GpB = lapply(UltraMammal, function(x) {
  compare.evol.rates(CarnMeanShapesTr$coords, x, gp.end, iter = 999, method = "simulation", print.progress = T)
})

Rate_Ratio<-unlist(lapply(Mammalcompareratesmulti.GpB, "[[", "sigma.d.ratio"))
Rate_P<-unlist(lapply(Mammalcompareratesmulti.GpB, "[[", "P.value"))
Rate_ES<-unlist(lapply(Mammalcompareratesmulti.GpB, "[[", "Z"))
Rate_Clade<-as.data.frame(matrix(unlist(lapply(Mammalcompareratesmulti.GpB, "[[", "sigma.d.gp")), nrow = 1000, ncol = 3, byrow = T))
Rate_Ppairwise<-as.data.frame(matrix(unlist(lapply(Mammalcompareratesmulti.GpB, "[[", "pairwise.pvalue")), nrow = 1000, ncol = 3, byrow = T))

#Collate
Rate_Data<-cbind.data.frame(Rate_Ratio,Rate_P,Rate_ES,Rate_Clade,Rate_Ppairwise)
colnames(Rate_Data)<-c("sigmadratio","P","Z","CanSigma", "DasSigma", "FelSigma", "CanDas_p", "CanFel_p", "DasFel_p")
write.csv(Rate_Data, "Rate_Data_CanFelDas-1K.csv")


##Compare rate metrics
setwd("/Users/home/Dropbox/Amherst/Project/Moduarity/MarsupialPlacental/3D/Paper/Submission/ResubResults")
#A
RateA<-read.csv("Rates/Rate_Data_CarDas-1K.csv", row.names = 1)
colMeans(RateA); apply(X = RateA, MARGIN = 2, FUN = median, na.rm = T)

apply(X = RateA, MARGIN = 2, FUN = std.error, na.rm = T)
apply(X = RateA, MARGIN = 2, FUN = quantile, probs=c(0.025,0.975), na.rm = T)

#B
RateB<-read.csv("Rates/Rate_Data_CanFelDas-1K.csv", row.names = 1)
colMeans(RateB); apply(X = RateB, MARGIN = 2, FUN = median, na.rm = T)

apply(X = RateB, MARGIN = 2, FUN = std.error, na.rm = T)
apply(X = RateB, MARGIN = 2, FUN = quantile, probs=c(0.025,0.975), na.rm = T)


###Plot Rates
library(ggplot2)
##Car-Das
Mamra<-read.csv("Rates/Rate_Data_CarDas-1Kstack.csv")
#Violin plots
ggplot(Mamra, aes(x=factor(Group), y=Rate, fill=factor(Group))) + geom_violin(width=1, trim=F) + labs(title="Rates", y = bquote('Brownian rate ('*~sigma^2*')'), x = "Clade", fill = "Clade") + theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "gray"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "gray"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_fill_discrete(name="Clade", breaks=c("Das", "Ecar"), labels=c("Dasyuromorphia", "Carnivora"))

##Das-Can-Fel
Mamra<-read.csv("Rates/Rate_Data_CanFelDas-1Kstack.csv")
#Violin plots
ggplot(Mamra, aes(x=factor(Group), y=Rate, fill=factor(Group))) + geom_violin(width=1, trim=F) + labs(title="Rates", y = bquote('Brownian rate ('*~sigma^2*')'), x = "Clade", fill = "Clade") + theme(panel.background = element_rect(fill = "white", colour = "black", size = 0.5, linetype = "solid"), panel.grid.major = element_line(size = 0.25, linetype = 'solid', colour = "gray"), panel.grid.minor = element_line(size = 0.25, linetype = 'solid',colour = "gray"), axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  scale_fill_discrete(name="Clade", breaks=c("Das", "Ecan", "Fel"), labels=c("Dasyuromorphia", "Caniformia", "Feliformia"))



####EMMLi Phylo####
library(EMMLi); library(paleotree); library(paleomorph)
setwd("/Users/home/Dropbox/Amherst/Project/Moduarity/MarsupialPlacental/3D/Paper/Submission/ResubResults/Writing/Final/GitHub/")
UltraMammal<-read.nexus(file = "UltraMammal1K_2019.nex")
CarnMeanData<-read.csv("Mean-Data_20210318.csv")
#Generate mean shapes above (line 1-51)

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

CarnMeanShapesTrEMMLi<-CarnMeanShapesTrEMMLi[,,as.character(KeepDas$PhyloID)]
CarnMeanShapesTrEMMLi<-CarnMeanShapesTrEMMLi[,,as.character(KeepCar$PhyloID)]
CarnMeanShapesTrEMMLi<-CarnMeanShapesTrEMMLi[,,as.character(KeepCan$PhyloID)]
CarnMeanShapesTrEMMLi<-CarnMeanShapesTrEMMLi[,,as.character(KeepFel$PhyloID)]

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


LMList<-read.csv("LMsmodel.csv")
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
save(Phylo.emmli.Das, file = "Ultrametric_PhyloemmliDas.rda")
save(Phylo.emmli.Car, file = "Ultrametric_PhyloemmliCar.rda")
save(Phylo.emmli.Can, file = "Ultrametric_PhyloemmliCan.rda")
save(Phylo.emmli.Fel, file = "Ultrametric_PhyloemmliFel.rda")


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
EffectQuantDas <- aaply(laply(ojEMMLiDas, as.matrix), c(2, 3), quantile, probs=c(0.025,0.975))
rownames(EffectMeanDas)<-rownames(Phylo.emmli.Das[[1]]$results)
write.csv(EffectMeanDas,"EMMLiDasResultsTableMean.csv")

#Carnivora
ojEMMLiCar<-list()
for(i in 1:1000){
  ojEMMLiCar[[i]]<-Phylo.emmli.Car[[i]]$results
  rownames(ojEMMLiCar[[i]]) <- c()
}

EffectMeanCar <- as.data.frame(aaply(laply(ojEMMLiCar, as.matrix), c(2, 3), mean))
EffectQuantCar <- aaply(laply(ojEMMLiCar, as.matrix), c(2, 3), quantile, probs=c(0.025,0.975))
rownames(EffectMeanCar)<-rownames(Phylo.emmli.Car[[1]]$results)
write.csv(EffectMeanCar,"EMMLiCarResultsTableMean.csv")

#Caniformia
ojEMMLiCan<-list()
for(i in 1:1000){
  ojEMMLiCan[[i]]<-Phylo.emmli.Can[[i]]$results
  rownames(ojEMMLiCan[[i]]) <- c()
}

EffectMeanCan <- as.data.frame(aaply(laply(ojEMMLiCan, as.matrix), c(2, 3), mean))
EffectQuantCan <- aaply(laply(ojEMMLiCan, as.matrix), c(2, 3), quantile, probs=c(0.025,0.975))
rownames(EffectMeanCan)<-rownames(Phylo.emmli.Can[[1]]$results)
write.csv(EffectMeanCan,"EMMLiCanResultsTableMean.csv")

#Feliformia
ojEMMLiFel<-list()
for(i in 1:1000){
  ojEMMLiFel[[i]]<-Phylo.emmli.Fel[[i]]$results
  rownames(ojEMMLiFel[[i]]) <- c()
}

EffectMeanFel <- as.data.frame(aaply(laply(ojEMMLiFel, as.matrix), c(2, 3), mean))
EffectQuantFel <- aaply(laply(ojEMMLiFel, as.matrix), c(2, 3), quantile, probs=c(0.025,0.975))
rownames(EffectMeanFel)<-rownames(Phylo.emmli.Fel[[1]]$results)
write.csv(EffectMeanFel,"EMMLiFelResultsTableMean.csv")


####EMMLi module consistency####
##Das
EMsmodulesDas<-list()
for(i in 1:1000){
  EMsmodulesDas[[i]]<-names(which(Phylo.emmli.Das[[i]]$results[,5]==0))
}
EMsmodulesDasEdit = lapply(EMsmodulesDas, function(x) {
  unique(gsub("^([^.]+.[^.]+).*","\\1",x))
})
table(unlist(EMsmodulesDasEdit))

##Car
EMsmodulesCar<-list()
for(i in 1:1000){
  EMsmodulesCar[[i]]<-names(which(Phylo.emmli.Car[[i]]$results[,5]==0))
}
EMsmodulesCarEdit = lapply(EMsmodulesCar, function(x) {
  unique(gsub("^([^.]+.[^.]+).*","\\1",x))
})
table(unlist(EMsmodulesCarEdit))

##Can
EMsmodulesCan<-list()
for(i in 1:1000){
  EMsmodulesCan[[i]]<-names(which(Phylo.emmli.Can[[i]]$results[,5]==0))
}
EMsmodulesCanEdit = lapply(EMsmodulesCan, function(x) {
  unique(gsub("^([^.]+.[^.]+).*","\\1",x))
})
table(unlist(EMsmodulesCanEdit))

##Fel
EMsmodulesFel<-list()
for(i in 1:1000){
  EMsmodulesFel[[i]]<-names(which(Phylo.emmli.Fel[[i]]$results[,5]==0))
}
EMsmodulesFelEdit = lapply(EMsmodulesFel, function(x) {
  unique(gsub("^([^.]+.[^.]+).*","\\1",x))
})
table(unlist(EMsmodulesFelEdit))


##Plotting
library(ggplot2)
#Collate the data
CanDasEm<-read.csv("/Users/home/Dropbox/Amherst/Project/Moduarity/MarsupialPlacental/3D/Paper/Submission/ResubResults/EMMLi/Mam_EMMLiSummaryPlot_wSize2.csv")

# Stacked + percent
ggplot(CanDasEm, aes(fill=Module, y=Frequency, x=Clade)) + 
  geom_bar(position="fill", stat="identity") +
  theme_bw()



####Zcr####
library(geiger); library(geomorph); library(paleomorph); library(EMMLi); library(phytools); library(plotrix); library(parallel)
setwd("/Users/home/Dropbox/Amherst/Project/Moduarity/MarsupialPlacental/3D/Paper/Submission/ResubResults/Writing/Final/GitHub/")

MandibleShape<-readland.tps("Individual-MarPlaJaw_20200805.tps", specID = "ID")
DataNak <-read.csv("Individual-Data_20200829_order.csv")
morphmean <- aggregate(CentroidSize~PhyloName,data=DataNak,FUN="mean",na.rm=TRUE,na.action=NULL)

Y.gpa<-gpagen(A = MandibleShape)

###Mean TPS file###
#Coordinates#
DataNak <-read.csv("Individual-Data_20200829_order.csv")
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


###Tree data engineering###
#Read distribution of trees
UltraMammal<-read.nexus(file = "UltraMammal1K_2019.nex")
CarnMeanData<-read.csv("Mean-Data_20210318.csv")
CarnMeanData$GroupA<-as.factor(CarnMeanData$GroupA)
CarnMeanData$GroupB<-as.factor(CarnMeanData$GroupB)

Y.gpa<-gpagen(A = CarnMeanShapes)
CarnMeanShapes<-Y.gpa$coords

#Retain taxa in tree
CarnMeanShapesTr<-CarnMeanShapes[,,dimnames(CarnMeanShapes)[[3]][dimnames(CarnMeanShapes)[[3]]%in%UltraMammal$tree_4444$tip.label]]

#Match TPS order to tree order
CarnMeanShapesTr<-CarnMeanShapesTr[,,(match(UltraMammal$tree_4444$tip.label, dimnames(CarnMeanShapesTr)[[3]]))]

#Procrustes and add back in mean centroid size
CarnMeanShapesTr<-gpagen(CarnMeanShapesTr)

LMList<-read.csv("LMsmodel.csv")
LMList[,1]<-as.character(LMList[,1])

#Change Model.8 to whichever model you want to test.
##Dasyuromorphia
CarnMeanShapesTrEMMLi<-CarnMeanShapesTr$coords

#Subset by clades
KeepDas<-CarnMeanData[CarnMeanData$GroupB=="Das",]

CarnMeanShapesTrEMMLiD<-CarnMeanShapesTrEMMLi[,,as.character(KeepDas$PhyloID)]

#Prune Trees
EMMLiDropTree<-lapply(UltraMammal,keep.tip,tip=dimnames(CarnMeanShapesTrEMMLiD)[[3]])
class(EMMLiDropTree)<-"multiPhylo"

mD.test = lapply(EMMLiDropTree, function(x) {
  phylo.modularity(A = CarnMeanShapesTrEMMLiD[1:11,,], partition.gp = LMList$Model.8[1:11],
                   phy = x, iter = 999, print.progress = T)
})



##Caniformes
CarnMeanShapesTrEMMLi<-CarnMeanShapesTr$coords

#Subset by clades
KeepCan<-CarnMeanData[CarnMeanData$GroupB=="Can",]

CarnMeanShapesTrEMMLiCa<-CarnMeanShapesTrEMMLi[,,as.character(KeepCan$PhyloID)]

#Prune Trees
EMMLiDropTree<-lapply(UltraMammal,keep.tip,tip=dimnames(CarnMeanShapesTrEMMLiCa)[[3]])
class(EMMLiDropTree)<-"multiPhylo"

mCa.test = lapply(EMMLiDropTree, function(x) {
  phylo.modularity(A = CarnMeanShapesTrEMMLiCa[1:11,,], partition.gp = LMList$Model.8[1:11],
                   phy = x, iter = 999, print.progress = T)
})


##Feliformes
CarnMeanShapesTrEMMLi<-CarnMeanShapesTr$coords

#Subset by clades
KeepFel<-CarnMeanData[CarnMeanData$GroupB=="Fel",]

CarnMeanShapesTrEMMLiF<-CarnMeanShapesTrEMMLi[,,as.character(KeepFel$PhyloID)]

#Prune Trees
EMMLiDropTree<-lapply(UltraMammal,keep.tip,tip=dimnames(CarnMeanShapesTrEMMLiF)[[3]])
class(EMMLiDropTree)<-"multiPhylo"

mF.test = lapply(EMMLiDropTree, function(x) {
  phylo.modularity(A = CarnMeanShapesTrEMMLiF[1:11,,], partition.gp = LMList$Model.8[1:11],
                   phy = x, iter = 999, print.progress = T)
})


##Carnivora
CarnMeanShapesTrEMMLi<-CarnMeanShapesTr$coords

#Subset by clades
KeepCar<-CarnMeanData[CarnMeanData$GroupA=="Car",]

CarnMeanShapesTrEMMLiC<-CarnMeanShapesTrEMMLi[,,as.character(KeepCar$PhyloID)]

#Prune Trees
EMMLiDropTree<-lapply(UltraMammal,keep.tip,tip=dimnames(CarnMeanShapesTrEMMLiC)[[3]])
class(EMMLiDropTree)<-"multiPhylo"

mC.test = lapply(EMMLiDropTree, function(x) {
  phylo.modularity(A = CarnMeanShapesTrEMMLiC[1:11,,], partition.gp = LMList$Model.8[1:11],
                   phy = x, iter = 999, print.progress = T)
})

#save(mD.test, file = "Das_8zCR_999.rda")
#save(mCa.test, file = "Can_8zCR_999.rda")
#save(mF.test, file = "Fel_8zCR_999.rda")
#save(mC.test, file = "Car_8zCR_999.rda")

#Combine and analyze
CRTest<-mapply(compare.CR, mD.test, mCa.test, mF.test, CR.null = F, SIMPLIFY=FALSE)

m8CompareSZ <- as.data.frame(matrix(unlist(lapply(CRTest,"[[", "sample.z")), nrow = 1000, ncol = 4, byrow = T))
m8ComparePZ <- as.data.frame(matrix(unlist(lapply(CRTest,"[[", "pairwise.z")), nrow = 1000, ncol = 16, byrow = T))[,c(2,3,4,7,8,12)]
m8ComparePP <- as.data.frame(matrix(unlist(lapply(CRTest,"[[", "pairwise.P")), nrow = 1000, ncol = 16, byrow = T))[,c(2,3,4,7,8,12)]

colMeans(m8ComparePP)
apply(X = m8ComparePP, MARGIN = 2, FUN = median, na.rm = T)

write.csv(m8CompareSZ,"m8CompareSZ.csv")
write.csv(m8ComparePZ,"m8ComparePZ.csv")
write.csv(m8ComparePP,"m8ComparePP.csv")

#CR outputs, change mD.test to other clades
mD.test.output <- as.data.frame(matrix(
  c(names(unlist(lapply(mD.test,"[[", "CR"))),
    unlist(lapply(mD.test,"[[", "CR")),
    unlist(lapply(mD.test,"[[", "P.value")),
    unlist(lapply(mD.test,"[[", "Z"))),
  nrow = 1000, ncol = 4, byrow = F))
write.csv(mD.test.output,"DasCRoutput_m8.csv") #Save

DasCR<-read.csv("zCR/m8/DasCRoutput_m8.csv",row.names = 1)
CanCR<-read.csv("zCR/m8/CanCRoutput_m8.csv",row.names = 1)
CarCR<-read.csv("zCR/m8/CarCRoutput_m8.csv",row.names = 1)
FelCR<-read.csv("zCR/m8/FelCRoutput_m8.csv",row.names = 1)

colMeans(CanCR);colMeans(CarCR);colMeans(DasCR);colMeans(FelCR)
apply(X = CanCR, MARGIN = 2, FUN = median, na.rm = T)
apply(X = CarCR, MARGIN = 2, FUN = median, na.rm = T)
apply(X = DasCR, MARGIN = 2, FUN = median, na.rm = T)
apply(X = FelCR, MARGIN = 2, FUN = median, na.rm = T)

apply(X = CanCR, MARGIN = 2, FUN = std.error, na.rm = T)
apply(X = CarCR, MARGIN = 2, FUN = std.error, na.rm = T)
apply(X = DasCR, MARGIN = 2, FUN = std.error, na.rm = T)
apply(X = FelCR, MARGIN = 2, FUN = std.error, na.rm = T)

PP<-read.csv("zCR/m8/m8ComparePP.csv",row.names = 1)
PZ<-read.csv("zCR/m8/m8ComparePZ.csv",row.names = 1)
SZ<-read.csv("zCR/m8/m8CompareSZ.csv",row.names = 1)


colMeans(PP);colMeans(PZ);colMeans(SZ)
apply(X = PP, MARGIN = 2, FUN = median, na.rm = T)
apply(X = PZ, MARGIN = 2, FUN = median, na.rm = T)
apply(X = SZ, MARGIN = 2, FUN = median, na.rm = T)

apply(X = PP, MARGIN = 2, FUN = std.error, na.rm = T)
apply(X = PZ, MARGIN = 2, FUN = std.error, na.rm = T)
apply(X = SZ, MARGIN = 2, FUN = std.error, na.rm = T)
