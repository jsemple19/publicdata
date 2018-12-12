library(gdata)
library(GenomicRanges)
library(BSgenome.Celegans.UCSC.ce11)
library(rtracklayer)
library(ggbio)


# get rex sites from Albritton paper
#Albritton_elife2017="https://elifesciences.org/download/aHR0cHM6Ly9jZG4uZWxpZmVzY2llbmNlcy5vcmcvYXJ0aWNsZXMvMjM2NDUvZWxpZmUtMjM2NDUtc3VwcDEtdjEuemlw/elife-23645-supp1-v1.zip?_hash=8w773VoPpu6mx1PUQppgaFNrEniZQDyGrnfBDR6QDt8%3D"
#download.file(Albritton_elife2017,destfile="./Albritton_elife2017.zip")
#unzip("./Albritton_elife2017.zip",files="Tables/SupplementaryFile1D.xls",junkpaths=T)
#rexSites<-read.xls("./SupplementaryFile1D.xls")

#read in the ahringer files' table
ahringerfiles<-read.table("./ahringer.txt", header=TRUE, stringsAsFactors = FALSE)
ahringerfiles
for (i in 1:dim(ahringerfiles)[1]) {
download.file(ahringerfiles$Sample[i],destfile=paste0(ahringerfiles$FileName[i],".bw"))
}

# change chr names to UCSC formats
rexSites$Chromosome<-gsub("CHROMOSOME_","chr",rexSites$Chromosome)

# make Granges
rexGR<-GRanges(seqnames=rexSites$Chromosome,ranges=IRanges(start=rexSites$Start,end=rexSites$End),
               strand="*", rexSites[,c("Rank","Previous.Name","strength.category")])

# get liftover chain file from UCSC
download.file("http://hgdownload.cse.ucsc.edu/goldenPath/ce10/liftOver/ce10ToCe11.over.chain.gz",
              destfile="./ce10ToCe11.over.chain.gz")
system("gunzip ./ce10ToCe11.over.chain.gz")
chainFile<-import.chain("./ce10ToCe11.over.chain")

#lift rex sites from ce10 to ce11
rexGR_ce11<-unlist(liftOver(rexGR,chain=chainFile))
# you need to make sure all chromosomes are present in the seqlevels and seqlengths (get from BSgenome)
seqlevels(rexGR_ce11)<-seqlevels(Celegans)
seqlengths(rexGR_ce11)<-seqlengths(Celegans)

rexGR_ce11$strength.category<-factor(rexGR_ce11$strength.category,levels=c("strong","intermediate","weak"))

# plot rex sites on chromosome layout
autoplot(rexGR_ce11,layout="karyogram",aes(color = strength.category,fill = strength.category)) +
  theme_alignment() +
  ggtitle("Strength of rex sites")