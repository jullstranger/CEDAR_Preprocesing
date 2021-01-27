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

tissue<-args[i + 1]
npcs<-args[i + 2]
npc<-as.numeric(args [i + 3])
pcs2remove<-args[i+4]
pheno.dir<-args[i+5]


pheno.file<-paste(pheno.dir, tissue , "_Log2_RSN_Batch_Corr_4Plink.txt" , sep="")
cov.file<-paste(pheno.dir, tissue , "_" , npcs, "_PCs_after_correction.txt", sep="")
pcs2remove<-as.numeric(unlist(strsplit(pcs2remove , ";" )))

pcs2use<-1:npc

if(length(pcs2remove) > 1 | (length(pcs2remove) == 1 &  pcs2remove != 0) ){
	print(paste("we will remove " , length(pcs2remove) , "PCs from the model" , sep=" "))
	pcs2use<-setdiff(pcs2use, pcs2remove)
}


pheno<-read.table(pheno.file , header = TRUE, check.names=FALSE)
cov<-read.table(cov.file , header = TRUE)

res<-pheno[, 1:2]

prev.prc<-0
done=0
to.do<-ncol(pheno) - 2
for( i in 3:ncol(pheno)){
	pcs <- paste("PC" , pcs2use , sep="")
	vars <- paste(pcs,  collapse="+")
	regression <- paste0("Pheno", " ~ ",  vars )
	
	model<-lm( as.formula(regression) , data = data.frame( Pheno = pheno[, i] , cov[, 3 : (2+npc)]) ,
			na.action = na.exclude)
	resi<- residuals(model)
	
	res<-cbind(res, resi)
	colnames(res)[ncol(res)] <- colnames(pheno)[i]
	
	done<-done+1
	prc<-as.integer(100 * done /to.do )
	
	if(prc %% 5 == 0 & prc != prev.prc){
		print(paste( prc, "% done", sep=" "))
	}
	prev.prc<-prc
}

file2write<-paste(pheno.dir,  tissue , "_Log2_RSN_Batch_Corr4_"  , npc, "PCs_4Plink.txt" , sep="")
write.table(file=file2write, res , row.names=FALSE, quote = FALSE)

###################################################
##### usable probes coordinates ###################



