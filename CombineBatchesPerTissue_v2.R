# TODO: Add comment
# 
# Author: julia
###############################################################################


args= commandArgs()
i = 1
while(args[i]!="--args"){
	cat(args[i])
	i=i+1
}
tissue<-args[i+1]

print(paste("Tissue" , tissue, "will be processed" , sep=" "))

library(lumi)
library(sva)


input.dirB2<-"/home/mass/GRD/joel/Data/GeneExpression/B2/"
input.dirB1<-"/home/mass/GRD/joel/Data/GeneExpression/B1/"
input.dirRS<-"/home/mass/GRD/joel/Data/GeneExpression/RedoneSamples/"

sentBar.file<-"/home/mass/GRD/joel/Data/QC_Analysis_B2_290816/Sample_SentBarc_DetGenes.txt"
pheno.file<-"/home/mass/GRD/joel/Data/GeneExpression/Cedar_Pheno_V3.txt" 
outlieFile<-"/home/mass/GRD/joel/Data/QC_Analysis_B2_290816/Outliers_B2_RS_B1.csv"
#probes.file<-"/home/mass/GRD/joel/Data/ReAnnotate/GoodProbesCoordsReannotated_X_Y.txt"
probes.file<-"/home/mass/GRD/joel/Data/ReAnnotate/GoodProbesCoordsReannotated_V2.txt"


grid.output.dir<- "/home/mass/GRD/joel/Research/GEN/UAG/CROHN/GE_Data_Log2_RSN_SBC/" 

if ( ! file.exists(grid.output.dir)){
	dir.create(grid.output.dir)
	print(paste( grid.output.dir , "is created" , sep=" " ))
}

grid.output.dir<-paste(grid.output.dir , tissue, "/" , sep="" )

if ( ! file.exists(grid.output.dir)){
	dir.create(grid.output.dir)
	print(paste( grid.output.dir , "is created" , sep=" " ))
}
output.dir<- "/home/mass/GRD/joel/Data/QC_Analysis_merged_170517/" 

if ( ! file.exists(output.dir)){
	dir.create(output.dir)
	print(paste( output.dir , "is created" , sep=" " ))
}


output.dir<-paste(output.dir , tissue, "/" , sep="")

if ( ! file.exists(output.dir)){
	dir.create(output.dir)
	print(paste(output.dir , "is created" , sep=" " ))
}


Outliers_Detector<-function(df){
	result<-rep(0,2)
	Q2<-quantile(df , na.rm=TRUE)[2]
	Q4<-quantile(df , na.rm=TRUE)[4]
	
	iqr<-IQR(df, na.rm=TRUE)
	
	result[1]<-Q2 - 1.5 * iqr
	result[2]<-Q4 + 1.5 * iqr
	
	return(result)
	
}


#reportsB2<-c( "FinalReport_V2301013_.txt"  , "FinalReport_051113_.txt" )
#reportRS<-"FinalReport_121213_.txt"
#
#
#report<-reportsB2[1]
#fileName<-paste(input.dirB2, report, sep="")
#assign( paste(report, "lumi" , sep="_") , lumiR.batch(fileName ))
#ge<-exprs(get(paste(report, "lumi" , sep="_")) )
#sd<-se.exprs(get(paste(report, "lumi" , sep="_")))
#det<-detection(get(paste(report, "lumi" , sep="_")))	
#contr<-controlData(get(paste(report, "lumi" , sep="_")))
#
#report<-reportsB2[2]
#fileName<-paste(input.dirB2, report, sep="")
#assign( paste(report, "lumi" , sep="_") , lumiR.batch(fileName ))
#ge<-cbind(ge, exprs(get(paste(report, "lumi" , sep="_")) ) )
#sd<-cbind(sd, se.exprs(get(paste(report, "lumi" , sep="_")) ) )
#det<-cbind(det, detection(get(paste(report, "lumi" , sep="_"))))
#contr<-cbind(contr, controlData(get(paste(report, "lumi" , sep="_")))[, -c(1,2)] )
#
#report<-reportRS
#fileName<-paste(input.dirRS, report , sep="")
#assign( paste(report, "lumi" , sep="_") , lumiR.batch(fileName ))
#ge<-cbind(ge, exprs(get(paste(report, "lumi" , sep="_")) ))
#sd<-cbind(sd, se.exprs(get(paste(report, "lumi" , sep="_")) ))
#det<-cbind(det, detection(get(paste(report, "lumi" , sep="_"))))
#contr<-cbind(contr, controlData(get(paste(report, "lumi" , sep="_")) ) [, -c(1,2)]  )
#
#
#ids<-grep( tissue , colnames(ge))
#ge<-ge[, ids]
#sd<-sd[, ids]
#det<-det[, ids]
#
#ids<-grep( tissue , colnames(contr))
#contr<-contr[ , c(1,2,ids)]
## check
#sum(colnames(ge) !=colnames(sd))
#sum(colnames(ge) !=colnames(det))
#sum(colnames(contr) [-c(1,2)] != colnames(ge))
#
#report<-paste("FinalReport", tissue, ".txt", sep="")
#fileName<-paste(input.dirB1, report, sep="")
#assign( paste(report, "lumi" , sep="_") , lumiR.batch(fileName ))
#ge<-cbind( exprs(get(paste(report, "lumi" , sep="_"))) ,  ge)
#sd<-cbind( se.exprs(get(paste(report, "lumi" , sep="_")) ) , sd)
#det<-cbind( detection(get(paste(report, "lumi" , sep="_"))) , det)
#contr<-cbind(controlData(get(paste(report, "lumi" , sep="_")))  , contr[, -c(1,2)])
#
#
#print( paste("there are " , ncol(ge) , "samples for tissue", tissue,   sep=" " ))
#
#sum(colnames(ge) !=colnames(sd))
#sum(colnames(ge) !=colnames(det))
#sum(colnames(contr) [-c(1,2)] != colnames(ge))

# it seems that algorithms doesn't work because I have same individuals included multiple time in my matrix
# I  should remove them before calculating the distances matrix and also before calculating 
# the principal component

grid.input.dir<-paste("/home/mass/GRD/joel/Data/Merged_Batches_291216/" , tissue, "/" , sep="" )

ge<-read.table(paste (grid.input.dir , "Matrix_GE_" , tissue , ".txt" , sep="" ) , header = TRUE , check.names = FALSE)
sd<-read.table(paste (grid.input.dir , "Matrix_SD_" , tissue , ".txt" , sep="" ) , header = TRUE , check.names = FALSE)
det<-read.table(paste (grid.input.dir , "Matrix_DET_" , tissue , ".txt" , sep="" ) , header = TRUE , check.names = FALSE)
contr<-read.table(paste (grid.input.dir , "Matrix_CONTR_" , tissue , ".txt" , sep="" ) , header = TRUE , check.names = FALSE)

##########################" thus before making the merged batch object, you should remove duplicated individuals

expressed<-function(v){
	good<-sum(v < 0.05 )
	return(good)
}

probes.expressed<-apply(det, 1, expressed )
bad<-which(probes.expressed < round( ncol(det)/4 ) )

expressed.per.ind<-colSums(det < 0.05)
## remove outlies and duplicates now, after this you can do transformation
getIPCnumber = function(IPC)
{
	template = NULL
	template = strsplit(IPC,split="PC")[[1]][2];
	if(length(grep("_" ,  template) ) == 1 ){
		template = strsplit(template,split="_")[[1]][1];
	}else{
		if(length(grep(pattern="bis|ter|quat|quin|sex",x=template)) == 1){
			template = strsplit(template,split="bis|ter|quat|quin|sex")[[1]][1];
		}
	}
	
	template = paste('IPC',template,sep="")
}

indis<-unlist(lapply (names(expressed.per.ind) , getIPCnumber ) )
expressed.per.ind.t<-data.frame(IID=indis, Name = names(expressed.per.ind),  ExprGenes = expressed.per.ind)


expressed.per.ind.t<-expressed.per.ind.t[order(expressed.per.ind.t$ExprGenes, decreasing=TRUE) , ]
duplis<-which(duplicated(expressed.per.ind.t$IID))

to.delete<-as.character(expressed.per.ind.t$Name[duplis ])
indxs<-which(colnames(ge) %in% to.delete )

ge<-ge[, -indxs]
sd<-sd[ , -indxs]
det<-det[ , -indxs]
indxs<-which(colnames(contr) %in% to.delete )
contr<-contr[ , -indxs]

print(paste ( "Original samples (" ,  tissue , " ) = " , length( indis ), " " , length(indxs) , " duplicated removed" , sep=" " ))

#############################################################################################################
merged.lumi<-new("LumiBatch", exprs = as.matrix(ge ) , se.exprs  =  as.matrix(sd ) ,
		detection = as.matrix(det ), controlData = contr  )

# this function doesn't work for CD4, could it be that I do have duplicated individualse here?
outlies<-detectOutlier(merged.lumi , ifPlot=FALSE)
oooo.sam<-names(outlies)[which(outlies==TRUE)]

indis<-colnames(ge)
ii<-substr(indis , 6, nchar(indis))

print(paste("Lumi has found" , length(oooo.sam) , "outliers" , sep=" "))

if(length(oooo.sam) > 0){
	bad.lumi<-which(colnames(ge) %in%  oooo.sam)
	cols<-rep( "olivedrab" , ncol(ge) )
	cols[bad.lumi]<-"maroon4"
	
	jpeg(paste(output.dir, "Samples_Relation_lumi_outl.jpeg" , sep=""), width=1024, height=1024)
	plotSampleRelation(as.matrix( ge ) , method= 'mds' , color=cols)
	dev.off()
}


outli.tab<-read.table(outlieFile , header = TRUE)
tis.outlies<-outli.tab$IID[grep (tissue , outli.tab$IID)]

common.outlies<-intersect(tis.outlies, oooo.sam)
print(paste("I have detected " , length(tis.outlies) , "outliers, common with lumi " ,  length(common.outlies) , sep=" "))
out.comb<-union(oooo.sam , outli.tab$IID[grep (tissue , outli.tab$IID)])
print(paste("Combined outliers " , length(out.comb)  , sep=" "))

bad.comb<-which(colnames(ge) %in%  out.comb)
cols<-rep( "olivedrab" , ncol(ge) )
cols[bad.comb]<-"maroon4"

jpeg(paste(output.dir, "Samples_Relation_comb_outl.jpeg" , sep=""), width=1024, height=1024)
plotSampleRelation(as.matrix( ge ), method= 'mds' , color=cols)
dev.off()

ge.p<-prcomp(t(ge) , center = TRUE, scale = TRUE)
p<-predict(ge.p)

# unfortunately coudn't use this library
# library(rgl)
# jpeg(paste(output.dir, "PCs_123.jpg" , sep=""), width=1024, height=1024)
# plot3d(p[, 1:3], col=rainbow(1:nrow(p)) )
# text3d(p[, 1:3],texts= ii)
# dev.off()
# check

sum(rownames(p) !=colnames(ge))

jpeg(paste(output.dir, "PCs_1_2.jpg" , sep=""), width=1024, height=1024)
plot(p[, 1] , p[,2] )
text(p[, 1] , p[,2] , ii , col=cols)
dev.off()

jpeg(paste(output.dir, "PCs_2_3.jpg" , sep=""), width=1024, height=1024)
plot(p[, 2] , p[,3] )
text(p[, 2] , p[,3] , ii, col= cols )
dev.off()

jpeg(paste(output.dir, "PCs_1_3.jpg" , sep=""), width=1024, height=1024)
plot(p[, 1] , p[,3] )
text(p[, 1] , p[,3] , ii , col=cols )
dev.off()

jpeg(paste(output.dir, "PCs_2_4.jpg" , sep=""), width=1024, height=1024)
plot(p[, 2] , p[,4] )
text(p[, 2] , p[,4] , ii , col=cols )
dev.off()


k.clust<-kmeans(t(ge) , centers=3)
dd<-data.frame( IID = names(k.clust$cluster) , Clust=k.clust$cluster)
sum(rownames(p) != dd$IID)

cols<-rep("chocolate1", nrow(p))
cols[dd[, 2] == 2]<- "maroon4"
cols[dd[, 2] == 3]<- "olivedrab"

jpeg(paste(output.dir, "PCs_1_3_k_means.jpg" , sep=""), width=1024, height=1024)
plot(p[, 1] , p[,3] )
text(p[, 1] , p[,3] , ii, col= cols )
dev.off()


#########################################################
probes.expressed<-apply(det, 1, expressed )

# yes, I should do it again, becaue I have deleted duplicates 
bad<-which(probes.expressed < round( ncol(det)/4 ) )
expressed.per.ind<-colSums(det < 0.05)
## remove outlies and duplicates now, after this you can do transformation
getIPCnumber = function(IPC)
{
	template = NULL
	template = strsplit(IPC,split="PC")[[1]][2];
	if(length(grep("_" ,  template) ) == 1 ){
		template = strsplit(template,split="_")[[1]][1];
	}else{
		if(length(grep(pattern="bis|ter|quat|quin|sex",x=template)) == 1){
			template = strsplit(template,split="bis|ter|quat|quin|sex")[[1]][1];
		}
	}
	
	template = paste('IPC',template,sep="")
}

indis<-unlist(lapply (names(expressed.per.ind) , getIPCnumber ) )
expressed.per.ind.t<-data.frame(IID=indis, Name = names(expressed.per.ind),  ExprGenes = expressed.per.ind)

expressed.per.ind.t<-expressed.per.ind.t[order(expressed.per.ind.t$ExprGenes, decreasing=TRUE) , ]
duplis<-which(duplicated(expressed.per.ind.t$IID))
# I shouldn't have duplicates here, no you shouldn't 

out.comb<-c(out.comb , "CD8IPC045_IPC046-A" , "CD8IPC045_IPC046-B" , "IL2IPC146")
to.delete<-intersect(expressed.per.ind.t$Name , out.comb)

##### outliers in per tissue 
print( paste ( "outliers (" , tissue, ") = " , length(to.delete) , sep=" ") )
# do background correction and log2 transform
ge.log2<-lumiT(merged.lumi, method = "log2")
ge.log2<-exprs(ge.log2)

indxs<-which(colnames(ge.log2) %in% to.delete )
ge.log2<-ge.log2[, - indxs]

indxs<-which(colnames(det) %in% to.delete )
det<-det[, - indxs]

# check
sum(colnames(ge.log2) != colnames(det))
ge.log2<-lumiN(ge.log2, method =  "rsn" )

indxs<-which(is.na(ge.log2) , arr.ind = TRUE)
# why should I do this , do you know?

if(nrow(indxs) > 0){
	ge.s<-ge.log2[-indxs[,1] , ]
}else{
	ge.s<-ge.log2
}

################################################################
barcodes<-read.table(sentBar.file , header  = TRUE)
dupli<-which(barcodes$Sen_Bar_Sec == "8986302039_D" )
if(length(dupli) > 0){
	barcodes<-barcodes[ - dupli , ]	
}

# now change the barcode for the sample 286 to 288, it will be correct (CD4 only)
tbl.m<-merge(barcodes ,expressed.per.ind.t , by.x="Sample_ID" , by.y="Name")

pheno.ss<-read.table(pheno.file , header = TRUE )
tbl.mm<-merge(tbl.m , pheno.ss , by.x="IID" , by.y="ind_IPC")

# extraction method stuff
extr.method.tis<-c("CD14","CD19", "IL","TR" , "RE")
if(tissue %in% extr.method.tis) {
	ipc<- as.integer(substr(tbl.mm$IID, 4 , 6) )
	extr.method<-rep(1, length(ipc))
	indxs.two<-which(ipc >=234)
	extr.method[indxs.two]<-2
	tbl.mm<-cbind(tbl.mm , ExtrMethod = extr.method)
}


# subset the phenotypes only for selected individuals
common.indis<-intersect(colnames(ge.s) , tbl.mm$Sample_ID)
tbl.mm<-tbl.mm[tbl.mm$Sample_ID %in% common.indis  , ]
ge.s<-ge.s[ , colnames(ge.s) %in% common.indis ]


rownames(tbl.mm)<-tbl.mm$Sample_ID
tbl.mm<-tbl.mm[ colnames(ge.s) , ]
# check
sum(colnames(ge.s) != tbl.mm$Sample_ID)

# "Sentrix_Barcode"  "    "Detected_Genes_0.05"    
# "Batch"               "ExprGenes"           "AGE"                
# "ind_ethnicity"       "ind_sex"             "BMI"                
# "ind_Smoker"          "Alcohol"             "ExtrMethod"  

effects<-c( "ExprGenes" ,"AGE" ,"Sentrix_Barcode", "ind_ethnicity"  , "ind_sex" , "ind_Smoker"  , "ExtrMethod" ) 

tbl.mm$Sentrix_Barcode<-as.factor(tbl.mm$Sentrix_Barcode)
tbl.mm$ind_ethnicity<-as.factor(tbl.mm$ind_ethnicity)
tbl.mm$ind_Smoker <- as.factor(tbl.mm$ind_Smoker)
tbl.mm$ind_sex <- as.factor(tbl.mm$ind_sex)
tbl.mm$Batch <- as.factor(tbl.mm$Batch)
if(sum(extr.method.tis == tissue) > 0 ){
	tbl.mm$ExtrMethod <- as.factor(tbl.mm$ExtrMethod)	
}


combat_edata<-ge.s	
indx<-which(colnames(tbl.mm) == "AGE")
age=tbl.mm[, indx]

indx<-which(colnames(tbl.mm) == "ExprGenes")
expr_g=tbl.mm[, indx]

indx<-which(colnames(tbl.mm) == "Batch")
batch=tbl.mm[, indx]


jpeg(paste(output.dir, "ExprsGenes_SntxBrc.jpg" , sep=""), width=1024, height=1024)
plot(tbl.mm[ , colnames(tbl.mm) %in%  c("Sentrix_Barcode" , "ExprGenes")])
dev.off()


jpeg(paste(output.dir, "ExprsGenes_Batch.jpg" , sep=""), width=1024, height=1024)
plot(tbl.mm[ , colnames(tbl.mm) %in%  c("Batch" , "ExprGenes")])
dev.off()

### now I want to do a trick, for each signleton array I find the brother having
### most similar number of expressed genes

per.batch<-aggregate(tbl.mm$Sentrix_Barcode, length , by=list(tbl.mm$Sentrix_Barcode) )
ii<-which(per.batch$x < 2)
b<-as.character(per.batch$Group.1[ii])
print( paste( "we have " , length(b) , "singletons") )
tbl.mm$Sentrix_Barcode<-as.character(tbl.mm$Sentrix_Barcode)
per.batch.expr<-aggregate(tbl.mm$ExprGenes, mean , by=list(tbl.mm$Sentrix_Barcode) )
iter<-0
while(length(b) > 0){
	s<-b[1]
	eg<-tbl.mm$ExprGenes[tbl.mm$Sentrix_Barcode == s]
	dist<-abs(eg - per.batch.expr$x)
	minBatch<-per.batch.expr$Group.1[which(dist == sort(dist)[2]) ]
	eg.n<-per.batch.expr$x[per.batch.expr$Group.1 == minBatch]
	
	print(paste("we want to assign array" , s , "eg" , eg , "to" , minBatch , "eg" , eg.n , sep=" " ))
	tbl.mm$Sentrix_Barcode[which(tbl.mm$Sentrix_Barcode == s)]<-minBatch
	per.batch<-aggregate(tbl.mm$Sentrix_Barcode, length , by=list(tbl.mm$Sentrix_Barcode) )
	ii<-which(per.batch$x < 2)
	b<-as.character(per.batch$Group.1[ii])
	
	tbl.mm$Sentrix_Barcode<-as.character(tbl.mm$Sentrix_Barcode)
	per.batch.expr<-aggregate(tbl.mm$ExprGenes, mean , by=list(tbl.mm$Sentrix_Barcode) )
	iter<-iter + 1
	print(paste("iter:" , iter ,  "we have another  " , length(b) , "singletons" , sep= " "))
	
}

tbl.mm$Sentrix_Barcode<-as.factor(tbl.mm$Sentrix_Barcode)
indx<-which(colnames(tbl.mm) == "Sentrix_Barcode")
sntrx=tbl.mm[, indx]

jpeg(paste(output.dir, "ExprsGenes_SntxBrc_After_SingletonRemove.jpg" , sep=""), width=1024, height=1024)
plot(tbl.mm[ , colnames(tbl.mm) %in%  c("Sentrix_Barcode" , "ExprGenes")])
dev.off()


modcombat = model.matrix( ~ AGE + ind_sex + ind_Smoker , data=tbl.mm)
combat_edata = ComBat(dat=combat_edata , batch = sntrx, mod=modcombat, par.prior=TRUE, prior.plot =FALSE )


Avg_Sign_AFC<-colMeans(combat_edata)
tbl.mm<-cbind(tbl.mm , Avg_Sign_AFC )

pdf(paste(output.dir, "Avg_Sign_SntxBrc_After_SingletonRemove_AFC.pdf" , sep=""), width=1024, height=1024)
plot(tbl.mm[ , colnames(tbl.mm) %in%  c("Sentrix_Barcode" , "Avg_Sign_AFC")])
dev.off()

jpeg(paste(output.dir, "Avg_Sign_SntxBrc_After_SingletonRemove_AFC.jpeg" , sep=""), width=1024, height=1024)
plot(tbl.mm[ , colnames(tbl.mm) %in%  c("Sentrix_Barcode" , "Avg_Sign_AFC")])
dev.off()

################################## why I do this? because I don't want imputed values in my dataset

excludedB2<-c( "FinalReport_V2301013_.txt_excluded_probesB2"  , "FinalReport_V2301013__excluded_probesB2.txt" )
excludedRS<-"FinalReport_121213__excluded_probesB2.txt"
excludedB1<-paste("FinalReport", tissue, "_excluded_probesB1.txt" , sep="")

exc.file<-paste(input.dirB2 , excludedB2[1] , sep="")
excluded<-read.table( exc.file, header = TRUE)
excluded<-excluded[excluded$Sample_ID %in% colnames(combat_edata) ,  ]
excluded<-excluded[excluded$Imputed==1 , ]
print(paste( "we have imputed" , nrow(excluded) , "probes-samples" , sep=" "))
exc.imp<-excluded

exc.file<-paste(input.dirB2 , excludedB2[2] , sep="")
excluded<-read.table( exc.file, header = TRUE)
excluded<-excluded[excluded$Sample_ID %in% colnames(combat_edata) ,  ]
excluded<-excluded[excluded$Imputed==1 , ]
print(paste( "we have imputed" , nrow(excluded) , "probes-samples" , sep=" "))
exc.imp<-rbind(exc.imp , excluded[ , -1])

exc.file<-paste(input.dirRS , excludedRS , sep="")
excluded<-read.table( exc.file, header = TRUE)
excluded<-excluded[excluded$Sample_ID %in% colnames(combat_edata) ,  ]
excluded<-excluded[excluded$Imputed==1 , ]
print(paste( "we have imputed" , nrow(excluded) , "probes-samples" , sep=" "))
exc.imp<-rbind(exc.imp , excluded[, -1])

exc.file<-paste(input.dirB1 , excludedB1 , sep="")
excluded<-read.table( exc.file, header = TRUE)
excluded<-excluded[excluded$Sample_ID %in% colnames(combat_edata) ,  ]
excluded<-excluded[excluded$Imputed==1 , ]
print(paste( "we have imputed" , nrow(excluded) , "probes-samples" , sep=" "))
exc.imp<-rbind(exc.imp , excluded [, -1])

duplis<-which(duplicated(exc.imp))
if(length(duplis) > 0){
	exc.imp <- exc.imp[-duplis , ]	
}

########################### Now create Principal component after exclusion of imputed probes
#### first I will exclude probes having imputed data in order to not influence the PCs

imp<-which(exc.imp$Imputed == 1)
probes2exclude<-unique(exc.imp$ProbeID[imp])

indx<-which(rownames(combat_edata) %in% probes2exclude)
m<-as.matrix(combat_edata[-indx ,  ])
##### I will try another method
ge.p<-prcomp(t(m) , center = TRUE, scale = TRUE)
p<-predict(ge.p)
#
#
jpeg(paste(output.dir, "PC_" , tissue, "_after_correction.jpeg" , sep=""))
plot(ge.p , main=paste("PCs for" , tissue , se=" ") )
dev.off()
#
jpeg(paste(output.dir, "PC_VarExpl" , tissue, "_after_correction.jpeg" , sep=""))
plot(ge.p , type="l" ,  main= paste(tissue , 
				" \n Cumulative Proportion PC10 = " , 
				round(summary(ge.p)$importance[3, 10], digits=2)  , sep="")  )
dev.off()
#
#
indis<-unlist(lapply (rownames(p) , getIPCnumber ) )
results<-data.frame(FID=indis, IID=indis, p)
npcs<-ncol(p)
file2write<-paste(grid.output.dir, tissue, "_", npcs, "_PCs_after_correction.txt" , sep="")
write.table(file = file2write, results , row.names = FALSE , quote = FALSE)

############################################################
# for absence of other ideas I will do the same for the dataset before batch correction
# exclude imputed probes here as well
m<-as.matrix(ge.s[-indx ,  ])
##### I will try another method
ge.p<-prcomp(t(m) , center = TRUE, scale = TRUE)
p<-predict(ge.p)
#
jpeg(paste(output.dir, "PC_" , tissue, "_not_corrected.jpeg" , sep=""))
plot(ge.p , main=paste("PCs for" , tissue , se=" ") )
dev.off()
#
jpeg(paste(output.dir, "PC_VarExpl" , tissue, "_not_corrected.jpeg" , sep=""))
plot(ge.p , type="l" ,  main= paste(tissue , 
				" \n Cumulative Proportion PC10 = " , 
				round(summary(ge.p)$importance[3, 10], digits=2)  , sep="")  )
dev.off()
indis<-unlist(lapply (rownames(p) , getIPCnumber ) )
results<-data.frame(FID=indis, IID=indis, p)
npcs<-ncol(p)

file2write<-paste(grid.output.dir, tissue, "_", npcs, "_PCs_not_corrected.txt" , sep="")
write.table(file = file2write, results , row.names = FALSE , quote = FALSE)

############################################################
###### no you should exclude probes which didn't pass the probes QC and you should also exclude probes 
# whcih didn't pass the mapping QC
# it should be the same can you check this , thus compare colnames combat with the colnames det

det<-det[, colnames (det) %in% colnames(combat_edata)]
probes.expressed<-apply(det, 1, expressed )
bad<-which(probes.expressed < round( ncol(det )/4 ) )


super.good<-read.table(probes.file , header = TRUE)

combat_edata<-combat_edata[ -bad , ]
ge.s<-ge.s[ -bad , ]


combat_edata.s<-combat_edata[rownames(combat_edata) %in% super.good$ProbeID,  ]
ge.s.s<-ge.s[rownames(ge.s) %in% super.good$ProbeID,  ]

sum(rownames(combat_edata.s) != rownames(ge.s.s))

indxs<-which(rownames(combat_edata.s) %in% probes2exclude)

for( indx in indxs){
	probe<-rownames(combat_edata.s)[indx]
	samples<-unique(as.character(exc.imp$Sample_ID[exc.imp$ProbeID == probe & exc.imp$Imputed == 1])) 
	print(paste( "for" , length(samples) , "samples NA is assigned for probe " , probe, sep=" " ) )
	
	row<-which(rownames(combat_edata.s) == probe)
	cols<-which(colnames(combat_edata.s) %in% samples )
	combat_edata.s[row, cols]<-NA
	ge.s.s[row, cols]<-NA
}

d <- t(combat_edata.s) 
probes<-colnames(d)
indis<-unlist(lapply (rownames(d) , getIPCnumber ) )

d<-data.frame(FID=indis ,  IID = indis ,  d)
colnames(d)<-c("FID" , "IID" , probes)

file2write<-paste(grid.output.dir, tissue, "_Log2_RSN_Batch_Corr_4Plink.txt" , sep="")
write.table(file=file2write, d , row.names=FALSE, quote = FALSE)

##################################################

d <- t(ge.s.s) 
d<-data.frame(FID=indis ,  IID = indis ,  d)
colnames(d)<-c("FID" , "IID" , probes)

file2write<-paste(grid.output.dir, tissue, "_Log2_RSN_Batch_Not_Corr_4Plink.txt" , sep="")
write.table(file=file2write, d , row.names=FALSE, quote = FALSE)


###########################################################################
# don't forget to eliminate principle components in another code

#### old and new comparison
#d <- t(combat_edata) 
#probes<-colnames(d)
#indis<-unlist(lapply (rownames(d) , getIPCnumber ) )
#
#d<-data.frame(FID=indis ,  IID = indis ,  d)
#colnames(d)<-c("FID" , "IID" , probes)
#
#d1<-read.table("~/Research/GEN/UAG/CROHN/GE_Data_Log2_RSN/CD14/CD14_Log2_RSN_Batch_Corr_4Plink.txt" , header = TRUE , check.names=FALSE)
#
#
## how do you want this compare?
## correlation for individuals and correlation between probes
#sum(colnames(d1) != colnames(d))
## 0
#sum(d1$FID != d$FID)
## 0
#probes.cor<-data.frame()
#prev.prc<-0
#done<-0
#nprobes<-ncol(d) - 2
#
#for(i in 3:ncol(d)){
#	corr<-cor(d[, i] , d1[, i] , use="pairwise.complete.obs")
#	probes.cor<-rbind(probes.cor , data.frame(colnames(d)[i] , corr ))
#	done<-done + 1
#	prc<-as.integer(done*100/nprobes)
#	if(prc %% 5 == 0 & prc !=prev.prc){
#		print(paste("done" , prc, "%", sep=" "))
#	}
#	prev.prc<-prc
#}
#
#quantile(probes.cor$corr)
##0%       25%       50%       75%      100% 
##0.6858768 0.9424909 0.9735009 0.9882948 0.9998151 
#
#
#ind.cor<-data.frame()
#prev.prc<-0
#done<-0
#nindis<-nrow(d) 
#
#for(i in 1:nrow(d)){
#	corr<-cor(d[i, 3:ncol(d)] , d1[i, 3:ncol(d)] , use="pairwise.complete.obs")
#	ind.cor<-rbind(ind.cor , data.frame(d$IID[i] , corr ))
#	done<-done + 1
#	prc<-as.integer(done*100/nindis)
#	if(prc %% 5 == 0 & prc !=prev.prc){
#		print(paste("done" , prc, "%", sep=" "))
#	}
#	prev.prc<-prc
#}

