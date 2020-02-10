library(bold)
library(tidyverse)
library(rentrez)

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

#Select rbcL
taxaBOLD_AquaNIS_rbcl_raw<-taxaBOLD_AquaNIS[taxaBOLD_AquaNIS$markercode=="rbcL"|taxaBOLD_AquaNIS$markercode=="rbcLa",]
taxaBOLD_EASIN_rbcl_raw<-taxaBOLD_EASIN[taxaBOLD_EASIN$markercode=="rbcL"|taxaBOLD_EASIN$markercode=="rbcLa",]

#Select rbcL>500
taxaBOLD_AquaNIS_rbcl<-taxaBOLD_AquaNIS_rbcl_raw[taxaBOLD_AquaNIS_rbcl_raw$number>500,]
taxaBOLD_EASIN_rbcl<-taxaBOLD_EASIN_rbcl_raw[taxaBOLD_EASIN_rbcl_raw$number>500,]

#Save and export dataframe downloaded from BOLD for rbcl>500
write.csv2(taxaBOLD_AquaNIS_rbcl, file="taxaBOLD_AquaNIS_rbcl.csv")
write.csv2(taxaBOLD_EASIN_rbcl, file="taxaBOLD_EASIN_rbcl.csv")


#Get stats for each species for rbcl
taxaBOLD_AquaNIS_rbcl_stats <- taxaBOLD_AquaNIS_rbcl %>%
  group_by(phylum_name, class_name, species_name, species_taxID)%>%
  summarize(Public_BOLD = table(species_name), GenBank_mining = sum(GenbankMined))

taxaBOLD_AquaNIS_rbcl_stats2<-bold_tax_id(id=taxaBOLD_AquaNIS_rbcl_stats$species_taxID, dataTypes = "stats")

taxaBOLD_AquaNIS_rbcl_statsfinal<-cbind.data.frame(taxaBOLD_AquaNIS_rbcl_stats, taxaBOLD_AquaNIS_rbcl_stats2)


taxaBOLD_EASIN_rbcl_stats <- taxaBOLD_EASIN_rbcl %>%
  group_by(phylum_name, class_name, species_name, species_taxID)%>%
  summarize(Public_BOLD = table(species_name), GenBank_mining = sum(GenbankMined))

taxaBOLD_EASIN_rbcl_stats2<-bold_tax_id(id=taxaBOLD_EASIN_rbcl_stats$species_taxID, dataTypes = "stats")

taxaBOLD_EASIN_rbcl_statsfinal<-cbind.data.frame(taxaBOLD_EASIN_rbcl_stats, taxaBOLD_EASIN_rbcl_stats2)

#Export excel
write.csv2(taxaBOLD_AquaNIS_rbcl_statsfinal, file="taxaBOLD_AquaNIS_stats_rbcl.csv")
write.csv2(taxaBOLD_EASIN_rbcl_statsfinal, file="taxaBOLD_EASIN_stats_rbcl.csv")


#Genbank Data Rbcl
#Rbcl, rubisco, according with wikipedia and genbank

#Download all genbank data rbcl>500
search_genbank_rbcl <- function(x){
  query2 <- paste(x, "[Organism] AND ((rbcl[Gene] OR rubisco[Gene]) AND (500[SLEN] : 3000[SLEN]))")
  rbcl<-entrez_search(db="nuccore", term=query2, retmax=10000)$count
}

GenBanktaxarbcl_AquaNIS<-lapply(taxa_AquaNIS$Taxa, search_genbank_rbcl)
GenBanktaxarbcl_AquaNISdf<-t(as.data.frame(GenBanktaxarbcl_AquaNIS))

write.csv2(GenBanktaxarbcl_AquaNISdf, file="taxaGenbank_AquaNIS_rbcl.csv")

GenBanktaxarbcl_EASIN<-lapply(taxa_EASIN$Taxa, search_genbank_rbcl)
GenBanktaxarbcl_EASINdf<-t(as.data.frame(GenBanktaxarbcl_EASIN))

write.csv2(GenBanktaxarbcl_EASINdf, file="taxaGenbank_EASIN_rbcl.csv")