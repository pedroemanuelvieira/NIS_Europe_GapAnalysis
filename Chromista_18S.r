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

#Select S18
taxaBOLD_AquaNIS_S18_raw<-taxaBOLD_AquaNIS[taxaBOLD_AquaNIS$markercode=="18S"|taxaBOLD_AquaNIS$markercode=="18Sa",]
taxaBOLD_EASIN_S18_raw<-taxaBOLD_EASIN[taxaBOLD_EASIN$markercode=="18S"|taxaBOLD_EASIN$markercode=="18Sa",]

#Select S18>500
taxaBOLD_AquaNIS_S18<-taxaBOLD_AquaNIS_S18_raw[taxaBOLD_AquaNIS_S18_raw$number>500,]
taxaBOLD_EASIN_S18<-taxaBOLD_EASIN_S18_raw[taxaBOLD_EASIN_S18_raw$number>500,]

#Save and export dataframe downloaded from BOLD for S18>500
write.csv2(taxaBOLD_AquaNIS_S18, file="taxaBOLD_AquaNIS_S18.csv")
write.csv2(taxaBOLD_EASIN_S18, file="taxaBOLD_EASIN_S18.csv")


#Get stats for each species for S18
taxaBOLD_AquaNIS_S18_stats <- taxaBOLD_AquaNIS_S18 %>%
  group_by(phylum_name, class_name, species_name, species_taxID)%>%
  summarize(Public_BOLD = table(species_name), GenBank_mining = sum(GenbankMined))

taxaBOLD_AquaNIS_S18_stats2<-bold_tax_id(id=taxaBOLD_AquaNIS_S18_stats$species_taxID, dataTypes = "stats")

taxaBOLD_AquaNIS_S18_statsfinal<-cbind.data.frame(taxaBOLD_AquaNIS_S18_stats, taxaBOLD_AquaNIS_S18_stats2)


taxaBOLD_EASIN_S18_stats <- taxaBOLD_EASIN_S18 %>%
  group_by(phylum_name, class_name, species_name, species_taxID)%>%
  summarize(Public_BOLD = table(species_name), GenBank_mining = sum(GenbankMined))

taxaBOLD_EASIN_S18_stats2<-bold_tax_id(id=taxaBOLD_EASIN_S18_stats$species_taxID, dataTypes = "stats")

taxaBOLD_EASIN_S18_statsfinal<-cbind.data.frame(taxaBOLD_EASIN_S18_stats, taxaBOLD_EASIN_S18_stats2)

#Export excel
write.csv2(taxaBOLD_AquaNIS_S18_statsfinal, file="taxaBOLD_AquaNIS_stats_18S.csv")
write.csv2(taxaBOLD_EASIN_S18_statsfinal, file="taxaBOLD_EASIN_stats_18S.csv")


#Genbank Data S18
#S18, rubisco, according with wikipedia and genbank

#Download all genbank data S18>500
search_genbank_S18 <- function(x){
  query2 <- paste(x, "[Organism] AND (18s ribosomal rna[Title] OR 18S rRNA[Title] OR 18S small subunit ribosomal RNA[Title] OR 18s ribosomal rna[Gene] OR 18S rRNA[Gene] OR 18S small subunit ribosomal RNA[Gene]) AND (500[SLEN]:5000000[SLEN]))")
  S18<-entrez_search(db="nuccore", term=query2, retmax=10000)$count
}

GenBanktaxaS18_AquaNIS<-lapply(taxa_AquaNIS$Taxa, search_genbank_S18)
GenBanktaxaS18_AquaNISdf<-t(as.data.frame(GenBanktaxaS18_AquaNIS))

write.csv2(GenBanktaxaS18_AquaNISdf, file="taxa_Genbank_AquaNIS_S18.csv")

GenBanktaxaS18_EASIN<-lapply(taxa_EASIN$Taxa, search_genbank_S18)
GenBanktaxaS18_EASINdf<-t(as.data.frame(GenBanktaxaS18_EASIN))

write.csv2(GenBanktaxaS18_EASINdf, file="taxa_Genbank_EASIN_S18.csv")
