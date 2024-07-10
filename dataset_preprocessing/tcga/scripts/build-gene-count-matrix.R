suppressPackageStartupMessages({
	library("readr")
	library("data.table")
	library("dplyr")
	library("tximport")
	library("GenomicFeatures")
})

file.not.empty<-function(filenames){file.info(filenames)$size != 0}

# IPUTS

args = commandArgs(trailingOnly=TRUE)

inpdir<-args[1]
outdir<-args[2]
metadat<-args[3]
cancer<-args[4]
txfile<-"/home/lmoral7/lmprojects/references/gencode_GRCh38.p13/gencode.v37.primary_assembly.annotation.tx2gene.tsv"

# MAIN 

print("Reading sample metadata...")
samp.info<-read.csv(metadat,header=T,row.names=1)
txg<-read.table(txfile,header=TRUE)


#Writing txi objects
################################################


print("=============Input directory============")
print(inpdir)

#Reading sample metadata

dirs<-list.dirs(inpdir,recursive=FALSE,full.names=F)
valid_ids<-lapply(dirs,function(d,...){
    file<-paste(inpdir,as.character(d),"quant.sf",sep="/")
    if(file.exists(file)&file.not.empty(file)){
        return(d)
    }
})%>%unlist()

print(paste("Reading",length(valid_ids),"sample files"))

stats<-data.frame(cancer=cancer,n_samples=nrow(samp.info),n_samples_salmon=length(valid_ids))
samp.info<-samp.info[valid_ids,]

# write metadata of valid files
metadat.out<-sub(".csv",".salmon.csv",metadat)
write.csv(samp.info,file=metadat.out)

files<-file.path(inpdir,rownames(samp.info),"quant.sf")

#Build tximport object
print("Reading files from input directory...")
txi<-tximport(files,type="salmon",tx2gene = txg, countsFromAbundance="no",dropInfReps=TRUE)
colnames(txi$counts)<-samp.info$submitter_id
colnames(txi$abundance)<-samp.info$submitter_id

cts<-as.data.frame(txi$counts)
tpm<-as.data.frame(txi$abundance)

#Save data
#######################################

print("Storing tximport object...")
objpath<-paste(outdir,"/",cancer,".gene.txiobject.RData",sep="")
save(txi,file=objpath)
print(paste("Writing object to:",objpath))

fpath<-paste(outdir,"/",cancer,".gene.counts.tsv",sep="")
write.table(cts,fpath)
print("Finished successfully!")

fpath<-paste(outdir,"/",cancer,".gene.tpm.tsv",sep="")
write.table(tpm,fpath)
print("Finished successfully!")