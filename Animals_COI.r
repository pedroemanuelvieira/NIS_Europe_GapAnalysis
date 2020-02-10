library(bold)
library(tidyverse)
library(rentrez)

###taxaCOI####
taxa_AquaNIS<-read.csv2("taxa_AquaNis.txt")
taxa_EASIN<-read.csv2("taxa_EASIN.txt")
#BOLD data

#Download bold data
taxaBOLD_AquaNIS1<-bold_seqspec(taxon=taxa_AquaNIS$Taxa[1:300], format = "tsv")
taxaBOLD_AquaNIS2<-bold_seqspec(taxon=taxa_AquaNIS$Taxa[301:600], format = "tsv")
taxaBOLD_AquaNIS3<-bold_seqspec(taxon=taxa_AquaNIS$Taxa[601:762], format = "tsv")
taxaBOLD_AquaNIS<-rbind(taxaBOLD_AquaNIS1, taxaBOLD_AquaNIS2,taxaBOLD_AquaNIS3)

taxaBOLD_EASIN1<-bold_seqspec(taxon=taxa_EASIN$Taxa[1:300], format = "tsv")
taxaBOLD_EASIN2<-bold_seqspec(taxon=taxa_EASIN$Taxa[301:600], format = "tsv")
taxaBOLD_EASIN3<-bold_seqspec(taxon=taxa_EASIN$Taxa[601:900], format = "tsv")
taxaBOLD_EASIN4<-bold_seqspec(taxon=taxa_EASIN$Taxa[901:1200], format = "tsv")
taxaBOLD_EASIN5<-bold_seqspec(taxon=taxa_EASIN$Taxa[1201:1219], format = "tsv")
taxaBOLD_EASIN<-rbind(taxaBOLD_EASIN1, taxaBOLD_EASIN2,taxaBOLD_EASIN3, taxaBOLD_EASIN4,taxaBOLD_EASIN5)

taxaBOLD_AquaNIS$number<-str_count(taxaBOLD_AquaNIS$nucleotides, pattern="[A-Z]")
taxaBOLD_EASIN$number<-str_count(taxaBOLD_EASIN$nucleotides, pattern="[A-Z]")

#Add information mined from genbank
taxaBOLD_AquaNIS$GenbankMined<-taxaBOLD_AquaNIS$institution_storing=="Mined from GenBank, NCBI"
taxaBOLD_EASIN$GenbankMined<-taxaBOLD_EASIN$institution_storing=="Mined from GenBank, NCBI"

#Select COI
taxaBOLD_AquaNIS_COI_raw<-taxaBOLD_AquaNIS[taxaBOLD_AquaNIS$markercode=="COI-5P",]
taxaBOLD_EASIN_COI_raw<-taxaBOLD_EASIN[taxaBOLD_EASIN$markercode=="COI-5P",]

#Select COI>500
taxaBOLD_AquaNIS_COI<-taxaBOLD_AquaNIS_COI_raw[taxaBOLD_AquaNIS_COI_raw$number>500,]
taxaBOLD_EASIN_COI<-taxaBOLD_EASIN_COI_raw[taxaBOLD_EASIN_COI_raw$number>500,]

#Save and export dataframe downloaded from BOLD for COI>500
write.csv2(taxaBOLD_AquaNIS_COI, file="taxaBOLD_AquaNIS_COI.csv")
write.csv2(taxaBOLD_EASIN_COI, file="taxaBOLD_EASIN_COI.csv")

#Get stats for each species for COI
taxaBOLD_AquaNIS_COI_stats <- taxaBOLD_AquaNIS_COI %>%
  group_by(phylum_name, class_name, species_name, species_taxID)%>%
  summarize(Public_BOLD = table(species_name), GenBank_mining = sum(GenbankMined))

taxaBOLD_AquaNIS_COI_stats2<-bold_tax_id(id=taxaBOLD_AquaNIS_COI_stats$species_taxID[ids], dataTypes = "stats")


#Export excel
write.csv2(taxaBOLD_AquaNIS_COI_stats, file="taxaBold_AquaNIS_COI_stats.csv")
write.csv2(taxaBOLD_AquaNIS_COI_stats2, file="taxaBold_AquaNIS_COI_stats2.csv")


taxaBOLD_AquaNIS_COI_statsfinal<-cbind.data.frame(taxaBOLD_AquaNIS_COI_stats, taxaBOLD_AquaNIS_COI_stats2)
#88
ids<-c(1:47,49:416, 418:489, 491:564)
taxaBOLD_EASIN_COI_stats <- taxaBOLD_EASIN_COI %>%
  group_by(phylum_name, class_name, species_name, species_taxID)%>%
  summarize(Public_BOLD = table(species_name), GenBank_mining = sum(GenbankMined))

taxaBOLD_EASIN_COI_stats2<-bold_tax_id(id=taxaBOLD_EASIN_COI_stats$species_taxID[ids], dataTypes = "stats")
#48,417,490
taxaBOLD_EASIN_COI_statsfinal<-cbind.data.frame(taxaBOLD_EASIN_COI_stats, taxaBOLD_EASIN_COI_stats2)

#Export excel
write.csv2(taxaBOLD_EASIN_COI_stats, file="taxaBold_EASIN_COI_stats.csv")
write.csv2(taxaBOLD_EASIN_COI_stats2, file="taxaBold_EASIN_COI_stats2.csv")


#Genbank Data

#Genbank Data COI>500
search_genbank_COI <- function(x){
  query2 <- paste(x, "[Organism] AND (((COI[Gene] OR CO1[Gene] OR COXI[Gene] OR COX1[Gene]) AND (500[SLEN]:3000[SLEN])) OR complete genome[All Fields] OR mitochondrial genome[All Fields])")
  COI<-entrez_search(db="nuccore", term=query2, retmax=10000)$count
}

GenBanktaxaCOI_AquaNIS<-lapply(taxa_AquaNIS$Taxa, search_genbank_COI)
GenBanktaxaCOI_AquaNISdf<-t(as.data.frame(GenBanktaxaCOI_AquaNIS))

write.csv2(GenBanktaxaCOI_AquaNISdf, file="taxaGenbank_AquaNIS_COI.csv")

GenBanktaxaCOI_EASIN<-lapply(taxa_EASIN$Taxa, search_genbank_COI)
GenBanktaxaCOI_EASINdf<-t(as.data.frame(GenBanktaxaCOI_EASIN))

write.csv2(GenBanktaxaCOI_EASINdf, file="taxaGenbank_EASIN_COI.csv")



