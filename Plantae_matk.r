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

#Select matk
taxaBOLD_AquaNIS_matk_raw<-taxaBOLD_AquaNIS[taxaBOLD_AquaNIS$markercode=="matK",]
taxaBOLD_EASIN_matk_raw<-taxaBOLD_EASIN[taxaBOLD_EASIN$markercode=="matK",]

#Select matk>500
taxaBOLD_AquaNIS_matk<-taxaBOLD_AquaNIS_matk_raw[taxaBOLD_AquaNIS_matk_raw$number>500,]
taxaBOLD_EASIN_matk<-taxaBOLD_EASIN_matk_raw[taxaBOLD_EASIN_matk_raw$number>500,]

#Save and export dataframe downloaded from BOLD for matk>500
write.csv2(taxaBOLD_AquaNIS_matk, file="taxaBOLD_AquaNIS_matk.csv")
write.csv2(taxaBOLD_EASIN_matk, file="taxaBOLD_EASIN_matk.csv")


#Get stats for each species for matk
taxaBOLD_AquaNIS_matk_stats <- taxaBOLD_AquaNIS_matk %>%
  group_by(phylum_name, class_name, species_name, species_taxID)%>%
  summarize(Public_BOLD = table(species_name), GenBank_mining = sum(GenbankMined))

taxaBOLD_AquaNIS_matk_stats2<-bold_tax_id(id=taxaBOLD_AquaNIS_matk_stats$species_taxID, dataTypes = "stats")

taxaBOLD_AquaNIS_matk_statsfinal<-cbind.data.frame(taxaBOLD_AquaNIS_matk_stats, taxaBOLD_AquaNIS_matk_stats2)


taxaBOLD_EASIN_matk_stats <- taxaBOLD_EASIN_matk %>%
  group_by(phylum_name, class_name, species_name, species_taxID)%>%
  summarize(Public_BOLD = table(species_name), GenBank_mining = sum(GenbankMined))

taxaBOLD_EASIN_matk_stats2<-bold_tax_id(id=taxaBOLD_EASIN_matk_stats$species_taxID, dataTypes = "stats")

taxaBOLD_EASIN_matk_statsfinal<-cbind.data.frame(taxaBOLD_EASIN_matk_stats, taxaBOLD_EASIN_matk_stats2)

#Export excel
write.csv2(taxaBOLD_AquaNIS_matk_statsfinal, file="taxaBOLD_AquaNIS_stats_matk.csv")
write.csv2(taxaBOLD_EASIN_matk_statsfinal, file="taxaBOLD_EASIN_stats_matk.csv")


#Genbank Data matk

#Download all genbank data matk>500
search_genbank_matk <- function(x){
  query2 <- paste(x, "[Organism] AND (matk[Gene] OR maturase K[Gene]) AND (500[SLEN] : 3000[SLEN]))")
  matk<-entrez_search(db="nuccore", term=query2, retmax=10000)$count
}

GenBanktaxamatk_AquaNIS<-lapply(taxa_AquaNIS$Taxa, search_genbank_matk)
GenBanktaxamatk_AquaNISdf<-t(as.data.frame(GenBanktaxamatk_AquaNIS))

write.csv2(GenBanktaxamatk_AquaNISdf, file="taxaGenbank_AquaNIS_matk.csv")

GenBanktaxamatk_EASIN<-lapply(taxa_EASIN$Taxa, search_genbank_matk)
GenBanktaxamatk_EASINdf<-t(as.data.frame(GenBanktaxamatk_EASIN))

write.csv2(GenBanktaxamatk_EASINdf, file="taxaGenbank_EASIN_matk.csv")
