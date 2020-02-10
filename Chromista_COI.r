library(bold)
library(tidyverse)
library(rentrez)

###taxaCOI####
taxa_AquaNIS<-read.csv2("taxa_AquaNis.txt")
taxa_EASIN<-read.csv2("taxa_EASIN.txt")
#BOLD data
#Download all bold data
taxaBOLD_AquaNIS<-bold_seqspec(taxon=taxa_AquaNIS$Taxa, format = "tsv")
taxaBOLD_EASIN<-bold_seqspec(taxon=taxa_EASIN$Taxa, format = "tsv")

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

taxaBOLD_AquaNIS_COI_stats2<-bold_tax_id(id=taxaBOLD_AquaNIS_COI_stats$species_taxID, dataTypes = "stats")

taxaBOLD_AquaNIS_COI_statsfinal<-cbind.data.frame(taxaBOLD_AquaNIS_COI_stats, taxaBOLD_AquaNIS_COI_stats2)


taxaBOLD_EASIN_COI_stats <- taxaBOLD_EASIN_COI %>%
  group_by(phylum_name, class_name, species_name, species_taxID)%>%
  summarize(Public_BOLD = table(species_name), GenBank_mining = sum(GenbankMined))

taxaBOLD_EASIN_COI_stats2<-bold_tax_id(id=taxaBOLD_EASIN_COI_stats$species_taxID, dataTypes = "stats")

taxaBOLD_EASIN_COI_statsfinal<-cbind.data.frame(taxaBOLD_EASIN_COI_stats, taxaBOLD_EASIN_COI_stats2)

#Export excel
write.csv2(taxaBOLD_AquaNIS_COI_statsfinal, file="taxaBOLD_AquaNIS_stats_COI.csv")
write.csv2(taxaBOLD_EASIN_COI_statsfinal, file="taxaBOLD_EASIN_stats_COI.csv")


#Genbank Data

#Genbank Data COI>500
search_genbank_COI <- function(x){
  query2 <- paste(x, "[Organism] AND (((COI[Gene] OR CO1[Gene] OR COXI[Gene] OR COX1[Gene]) AND (500[SLEN]:3000[SLEN])) OR 'complete genome'[All Fields] OR 'mitochondrial genome'[All Fields])")
  COI<-entrez_search(db="nuccore", term=query2, retmax=10000)$count
}

GenBanktaxaCOI_AquaNIS<-lapply(taxa_AquaNIS$Taxa, search_genbank_COI)
GenBanktaxaCOI_AquaNISdf<-t(as.data.frame(GenBanktaxaCOI_AquaNIS))

write.csv2(GenBanktaxaCOI_AquaNISdf, file="taxaGenbank_AquaNIS_COI.csv")

GenBanktaxaCOI_EASIN<-lapply(taxa_EASIN$Taxa, search_genbank_COI)
GenBanktaxaCOI_EASINdf<-t(as.data.frame(GenBanktaxaCOI_EASIN))

write.csv2(GenBanktaxaCOI_EASINdf, file="taxaGenbank_EASIN_COI.csv")



