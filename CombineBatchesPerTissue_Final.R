# TODO: Add comment
# 
# Author: julia
###############################################################################
##
## can you make two datasets one with X and Y chromosomes and other with all of them

#TISSUE=$1
#SCRATCH=/scratch/ulg/genan/joel/

#GE_DIR=~/Data/GeneExpression_DB/
#INPUT_DIR=~/Data/Merged_Batches_291216/${TISSUE}/
#GRID_OUT_DIR=${SCRATCH}GE_Data_Log2_RSN_SBC_v2/
#OUTPUT_DIR=~/Data/GE_Analysis_260717/
		
		
args= commandArgs()
i = 1
while(args[i]!="--args"){
	cat(args[i])
	i=i+1
}
tissue<-args[ i + 1]
ge.dir<-args[ i + 2]
input.dir<-args[ i + 3]
grid.output.dir<-args[ i + 4 ]
output.dir<-args[ i + 5]

print(paste("Tissue" , tissue, "will be processed" , sep=" "))

library(lumi)
library(sva)


input.dirB2<- paste(ge.dir , "B2/" , sep="")
input.dirB1<- paste(ge.dir , "B1/" , sep="")
input.dirRS<- paste (ge.dir, "RedoneSamples/" , sep="")

sentBar.file<-paste( ge.dir, "Sample_SentBarc_DetGenes.txt" , sep="")
pheno.file<- paste( ge.dir, "Cedar_Pheno_V3.txt" , sep="") 
outlieFile<- paste( ge.dir, "Outliers_B2_RS_B1.csv" , sep="")

cluster.outlie.file<-paste(ge.dir , "Cluster_Outliers.txt" , sep="")
swapped.samples.file<-paste(ge.dir , "swapped_samples.txt" , sep="")

probes.x.y.file<-"~/Data/ReAnnotate/GoodProbesCoordsReannotated_X_Y.txt"
probes.file<-"~/Data/ReAnnotate/GoodProbesCoordsReannotated_V2.txt"




if ( ! file.exists(grid.output.dir)){
	dir.create(grid.output.dir)
	print(paste( grid.output.dir , "is created" , sep=" " ))
}

grid.output.dir<-paste(grid.output.dir , tissue, "/" , sep="" )

if ( ! file.exists(grid.output.dir)){
	dir.create(grid.output.dir)
	print(paste( grid.output.dir , "is created" , sep=" " ))
}


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



ge<-read.table(paste (input.dir , "Matrix_GE_" , tissue , ".txt" , sep="" ) , header = TRUE , check.names = FALSE)
sd<-read.table(paste (input.dir , "Matrix_SD_" , tissue , ".txt" , sep="" ) , header = TRUE , check.names = FALSE)
det<-read.table(paste (input.dir , "Matrix_DET_" , tissue , ".txt" , sep="" ) , header = TRUE , check.names = FALSE)
contr<-read.table(paste (input.dir , "Matrix_CONTR_" , tissue , ".txt" , sep="" ) , header = TRUE , check.names = FALSE)

##########################" thus before making the merged batch object, you should remove duplicated individuals

expressed<-function(v){
	# detected with p-value < 0.05
	good<-sum(v < 0.05 )
	return(good)
}


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
if(tissue == "PLA"){
	to.delete<-c(to.delete , "PLAIPC427")
}
indxs<-which(colnames(ge) %in% to.delete )

if(length(indxs) > 0){
	ge<-ge[, -indxs]
	sd<-sd[ , -indxs]
	det<-det[ , -indxs]
	indxs<-which(colnames(contr) %in% to.delete )
	contr<-contr[ , -indxs]
	print(paste ( "Original sample number (" ,  tissue , " ) = " , 
					length( indis ), " " , length(indxs) , " duplicated removed" , sep=" " ))	
} 



#############################################################################################################
### we are merging data back without duplicates, should I remove simblings from here as well ?
### no, I will not , I will keep them untile the latest dataset is produced

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
tis.outlies<-as.character( outli.tab$IID[grep (tissue , outli.tab$IID)] )

common.outlies<-intersect(tis.outlies, oooo.sam)

print(paste("I have detected " , length(tis.outlies) , "outliers, common with lumi " ,  length(common.outlies) , sep=" "))
out.comb<-union(oooo.sam , tis.outlies )
print(paste("Combined outliers " , length(out.comb)  , sep=" "))

bad.comb<-which(colnames(ge) %in%  out.comb)
cols<-rep( "olivedrab" , ncol(ge) )
cols[bad.comb]<-"maroon4"

jpeg(paste(output.dir, "Samples_Relation_comb_outl.jpeg" , sep=""), width=1024, height=1024)
plotSampleRelation(as.matrix( ge ), method= 'mds' , color=cols)
dev.off()

ge.p<-prcomp(t(ge) , center = TRUE, scale = TRUE)
p<-predict(ge.p)



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


k.clust<-kmeans(t(ge) , centers=2)
dd<-data.frame( IID = names(k.clust$cluster) , Clust=k.clust$cluster)
sum(rownames(p) != dd$IID)

cols<-rep("chocolate1", nrow(p))
cols[dd[, 2] == 2]<- "maroon4"


jpeg(paste(output.dir, "PCs_1_3_k_means.jpg" , sep=""), width=1024, height=1024)
plot(p[, 1] , p[,3] )
text(p[, 1] , p[,3] , ii, col= cols )
dev.off()


#########################################################
# yes, I should do it again, becaue I have deleted duplicates 

expressed.per.ind<-colSums(det < 0.05)
## remove outlies and duplicates now, after this you can do transformation

indis<-unlist(lapply (names(expressed.per.ind) , getIPCnumber ) )
expressed.per.ind.t<-data.frame(IID=indis, Name = names(expressed.per.ind),  ExprGenes = expressed.per.ind)

expressed.per.ind.t<-expressed.per.ind.t[order(expressed.per.ind.t$ExprGenes, decreasing=TRUE) , ]
duplis<-which(duplicated(expressed.per.ind.t$IID))
# I shouldn't have duplicates here, no you shouldn't 
print(paste("there should be no duplicates and we have " , length(duplis) , "duplicates" , sep=" " ))

swapped.samples<-read.table(swapped.samples.file , header = TRUE , stringsAsFactors = FALSE)
cluster.outliers<-read.table(cluster.outlie.file, header = TRUE , stringsAsFactors = FALSE)
wierd.samples<-c("CD8IPC045_IPC046-A" , "CD8IPC045_IPC046-B" , "IL2IPC146")

out.comb<-c(out.comb , wierd.samples ,  swapped.samples$S1 , swapped.samples$S2 , cluster.outliers$x)

to.delete<-intersect(expressed.per.ind.t$Name , out.comb)
keept<-setdiff(expressed.per.ind.t$Name , out.comb)

##### outliers in per tissue 
print( paste ( "outliers (" , tissue, ") = " , length(to.delete) , "keept" , length(keept) ,  sep=" ") )
##### I  have to create new merged lumi package after removing of outliers
# do background correction and log2 transform

indxs<-which(colnames(ge) %in% to.delete )

if(length(indxs) > 0){
	print(paste ( "Previous sample number (" ,  tissue , " ) = " , 
					length( indis ), " " , length(indxs) , " outlies will be removed" , sep=" " ))	
	ge<-ge[, -indxs]
	sd<-sd[ , -indxs]
	det<-det[ , -indxs]
	indxs<-which(colnames(contr) %in% to.delete )
	contr<-contr[ , -indxs]
	merged.lumi<-new("LumiBatch", exprs = as.matrix(ge ) , se.exprs  =  as.matrix(sd ) ,
			detection = as.matrix(det ), controlData = contr  )
} 


# do log2 transform
ge.log2<-lumiT(merged.lumi, method = "log2")
ge.log2<-exprs(ge.log2)

# check
print( paste(" should be zero" , sum(colnames(ge.log2) != colnames(det)) , sep=" ")  )
# do robust spline normalization
ge.log2<-lumiN(ge.log2, method =  "rsn" )

indxs<-which(is.na(ge.log2) , arr.ind = TRUE)
# why should I do this , do you know?
ge.s<-ge.log2

################################################################
barcodes<-read.table(sentBar.file , header  = TRUE , stringsAsFactors = FALSE)
dupli<-which(barcodes$Sen_Bar_Sec == "8986302039_D" )
if(length(dupli) > 0){
	barcodes<-barcodes[ - dupli , ]	
}

# now change the barcode for the sample 286 to 288, it will be correct (CD4 only)
tbl.m<-merge(barcodes , expressed.per.ind.t , by.x="Sample_ID" , by.y="Name")

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
print(paste("should be zero" , sum(colnames(ge.s) != tbl.mm$Sample_ID)  , sep=" " ) ) 

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
if(length(indx) > 0){
	m<-as.matrix(combat_edata[-indx ,  ])	
	print(paste("we have excluded " , length(indx) , "imputed probes from PC calculation" , sep=" "))
}else{
	m<-as.matrix(combat_edata)
}

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
if(length(indx) > 0 ){
	m<-as.matrix(ge.s[-indx ,  ])	
}else{
	m<-as.matrix(ge.s)
}

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
combat_edata<-combat_edata[ -bad , ]
ge.s<-ge.s[ -bad , ]

super.good<-read.table(probes.file , header = TRUE , stringsAsFactors = FALSE)

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

super.good.x.y<-read.table(probes.x.y.file , header = TRUE, stringsAsFactors = FALSE)
super.good<-rbind(super.good , super.good.x.y)

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

file2write<-paste(grid.output.dir, tissue, "_Log2_RSN_Batch_Corr_4Plink_X_Y_Incl.txt" , sep="")
write.table(file=file2write, d , row.names=FALSE, quote = FALSE)

##################################################

d <- t(ge.s.s) 
d<-data.frame(FID=indis ,  IID = indis ,  d)
colnames(d)<-c("FID" , "IID" , probes)

file2write<-paste(grid.output.dir, tissue, "_Log2_RSN_Batch_Not_Corr_4Plink_X_Y_Incl.txt" , sep="")
write.table(file=file2write, d , row.names=FALSE, quote = FALSE)


###########################################################################

# don't forget to eliminate principle components in another code




