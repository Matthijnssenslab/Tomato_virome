library(tidyverse)
library(ggplot2)
library(tibble)
library(plyr)
library(dplyr)
library(reshape)
library(reshape2)
library(readxl)
library(devtools)
library(taxonomizr)
library(stringr)
library(openxlsx)
library(RColorBrewer)
library(grid)
library(readr)
library(data.table)
library(Polychrome)

### BUILD NR TAXONOMY OVERVIEW (One-off)###
##Before starting the R script run:
##1. ktClassifyBLAST NR_contigs_500_95-85.m8 -o NR_contigs_500_95-85.diamond.tab
##2. sort -k1,1 -k12,12gr -k11,11g -k3,3gr NR_contigs_500_95-85.m8 > test1.txt 
##   sort -u -k1,1 test1.txt > NR_contigs_500_95-85_diamond-best-hit.txt

# Update taxonomy database
getAccession2taxid(baseUrl='https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/')
prepareDatabase('accessionTaxa.sql',url= 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz')

#Import the sorted NR m8 file and name the columns
NR <- read_delim("NR_contigs_500_95-85_diamond-best-hit.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, col_types = cols())
colnames(NR) <- c('qseqid', 'sseqid', 'pident', 'align_length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')

#Add contig length to the overview table
words <- as.numeric(word(NR$qseqid, 4, sep = "_"))
length <- data.frame(matrix(unlist(words), nrow=length(words), byrow=T))
colnames(length) = c("Contig_length")
NR <- cbind(NR, length)

#Import the diamond.tab file and name the columns
taxdf <- read_delim("NR_contigs_500_95-85.diamond.tab", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, skip = 1, col_types = cols())
colnames(taxdf) <- c('qseqid', 'taxID', 'Avg. log e-value')

#Merge the dataframes by contig name
new <- merge(NR, taxdf, by = "qseqid", all.x = TRUE, all.y = TRUE)

#Attach taxonomy for the NR contigs
df <- as.data.frame(getTaxonomy(ids = new$taxID, sqlFile = "accessionTaxa.sql", desiredTaxa = c('superkingdom', 'phylum', 'class', 'order', 'family', 'subfamily', 'genus', 'species')), stringsAsFactors = FALSE)
overview <- cbind(new, df)

#Export the overview table
write.xlsx(overview, "NR_contigs_Diamond_overview.xlsx")


###ORGANISE DATA###
##Import the overview table created before

NR_contigs_Diamond_overview <- read_excel("NR_contigs_Diamond_overview.xlsx")

#Perform the following steps once to replace empty cells with "Unclassified". Repeat only if NR list is updated.
NR_contigs_Diamond_overview$superkingdom <- as.character(NR_contigs_Diamond_overview$superkingdom)
NR_contigs_Diamond_overview$phylum <- as.character(NR_contigs_Diamond_overview$phylum)
NR_contigs_Diamond_overview$class <- as.character(NR_contigs_Diamond_overview$class)
NR_contigs_Diamond_overview$order <- as.character(NR_contigs_Diamond_overview$order)
NR_contigs_Diamond_overview$family <- as.character(NR_contigs_Diamond_overview$family)
NR_contigs_Diamond_overview$subfamily <- as.character(NR_contigs_Diamond_overview$subfamily)
NR_contigs_Diamond_overview$genus <- as.character(NR_contigs_Diamond_overview$genus)
NR_contigs_Diamond_overview$species <- as.character(NR_contigs_Diamond_overview$species)

NR_contigs_Diamond_overview$superkingdom[is.na(NR_contigs_Diamond_overview$superkingdom)] <- "Unclassified"
NR_contigs_Diamond_overview$phylum[is.na(NR_contigs_Diamond_overview$phylum)] <- "Unclassified"
NR_contigs_Diamond_overview$class[is.na(NR_contigs_Diamond_overview$class)] <- "Unclassified"
NR_contigs_Diamond_overview$order[is.na(NR_contigs_Diamond_overview$order)] <- "Unclassified"
NR_contigs_Diamond_overview$family[is.na(NR_contigs_Diamond_overview$family)] <- "Unclassified"
NR_contigs_Diamond_overview$subfamily[is.na(NR_contigs_Diamond_overview$subfamily)] <- "Unclassified"
NR_contigs_Diamond_overview$genus[is.na(NR_contigs_Diamond_overview$genus)] <- "Unclassified"
NR_contigs_Diamond_overview$species[is.na(NR_contigs_Diamond_overview$species)] <- "Unclassified"

write.xlsx(NR_contigs_Diamond_overview, "NR_contigs_Diamond_overview.xlsx")

##Load genomad results

#Before moving to the next command, edit the _virus_summary.tsv file so that you have separate columns for each taxonomy step. 
#Name the first column qseqid and name the taxonomy steps as indicated below (eg gen_superkingdom).
genomad <- read_delim("NR_contigs_500_95-85_virus_summary_edit.txt", "\t", escape_double = FALSE, col_names = TRUE, col_types = cols())

#Replace empty cells with Unclassified, add genus and species columns with Unclassified
genomad$gen_superkingdom <- as.character(genomad$gen_superkingdom)
genomad$gen_phylum <- as.character(genomad$gen_phylum)
genomad$gen_class <- as.character(genomad$gen_class)
genomad$gen_order <- as.character(genomad$gen_order)
genomad$gen_family <- as.character(genomad$gen_family)
genomad$gen_superkingdom[is.na(genomad$gen_superkingdom)] <- "Unclassified"
genomad$gen_phylum[is.na(genomad$gen_phylum)] <- "Unclassified"
genomad$gen_class[is.na(genomad$gen_class)] <- "Unclassified"
genomad$gen_order[is.na(genomad$gen_order)] <- "Unclassified"
genomad$gen_family[is.na(genomad$gen_family)] <- "Unclassified"

genomad <- mutate(genomad, gen_genus = "Unclassified")
genomad <- mutate(genomad, gen_species = "Unclassified")

#Merge genomad and Diamond results
NR_contigs_DiaGen <- merge(NR_contigs_Diamond_overview, genomad, by = "qseqid", all = TRUE)

##Create subsets to merge Genomad and Diamond taxonomy appropriately
#Genomad overrides Diamond taxonomy for misclassified prokaryotic viruses. In the two datasets that translates to: 
#a) contigs classified by Diamond as: Unclassified and Bacteria, 
#b) contigs classified by Genomad as: Caudoviricetes (host: bacteria)
caudo <- filter(NR_contigs_DiaGen, gen_class == 'Caudoviricetes')

Viral_DiaGen <- filter(NR_contigs_DiaGen, gen_superkingdom == 'Viruses')
bacGen <- filter(Viral_DiaGen, superkingdom == 'Bacteria')
duprows <- bacGen$qseqid %in% caudo$qseqid
newtest <- rbind(caudo, bacGen[!duprows,])

uncGen <- filter(Viral_DiaGen, superkingdom == 'Unclassified')
duprows <- uncGen$qseqid %in% newtest$qseqid
newtest <- rbind(newtest, uncGen[!duprows,])

#Create row with Override information
NR_contigs_DiaGen <- mutate(NR_contigs_DiaGen, GenOverride = ifelse(NR_contigs_DiaGen$qseqid %in% newtest$qseqid , 'yes', 'no'))

#Create final NR table where Genomad taxonomy overrides Diamond for misclassified prokaryotic viruses
NR_contigs_DiaGen <- mutate(NR_contigs_DiaGen, Final_Superkingdom = ifelse(NR_contigs_DiaGen$GenOverride == 'yes', NR_contigs_DiaGen$gen_superkingdom, NR_contigs_DiaGen$superkingdom))
NR_contigs_DiaGen <- mutate(NR_contigs_DiaGen, Final_Phylum = ifelse(NR_contigs_DiaGen$GenOverride == 'yes', NR_contigs_DiaGen$gen_phylum, NR_contigs_DiaGen$phylum))
NR_contigs_DiaGen <- mutate(NR_contigs_DiaGen, Final_Class = ifelse(NR_contigs_DiaGen$GenOverride == 'yes', NR_contigs_DiaGen$gen_class, NR_contigs_DiaGen$class))
NR_contigs_DiaGen <- mutate(NR_contigs_DiaGen, Final_Order = ifelse(NR_contigs_DiaGen$GenOverride == 'yes', NR_contigs_DiaGen$gen_order, NR_contigs_DiaGen$order))
NR_contigs_DiaGen <- mutate(NR_contigs_DiaGen, Final_Family = ifelse(NR_contigs_DiaGen$GenOverride == 'yes', NR_contigs_DiaGen$gen_family, NR_contigs_DiaGen$family))
NR_contigs_DiaGen <- mutate(NR_contigs_DiaGen, Final_Genus = ifelse(NR_contigs_DiaGen$GenOverride == 'yes', NR_contigs_DiaGen$gen_genus, NR_contigs_DiaGen$genus))
NR_contigs_DiaGen <- mutate(NR_contigs_DiaGen, Final_Species = ifelse(NR_contigs_DiaGen$GenOverride == 'yes', NR_contigs_DiaGen$gen_species, NR_contigs_DiaGen$species))

NR_contigs_Final <- subset(NR_contigs_DiaGen, select = -c(16, 17, 18, 19, 20, 21, 22, 23, 33, 34, 35, 36, 37, 38, 39, 40, 41))

##Our focus is eukaryotic viruses. ICTV and ViralZone was used to characterize viral findings as Eukaryotic or Prokaryotic viruses based on the hosts of the identified taxa. Genome composition information were also kept.  
#Load the table with the family-host association for the viral findings
FamilyHost <- read_excel("family-host_association.xlsx")

##Create overall findings column
#Make vector for eukaryotic viruses
Euk_FamilyHost <- filter(FamilyHost, virus == "Eukaryotic")
euk_family <- c(Euk_FamilyHost$family)

#Create the overall findings column 
NR_contigs_Final <- mutate(NR_contigs_Final, Findings = ifelse(NR_contigs_Final$Final_Superkingdom == 'Eukaryota', 'Eukaryota', ifelse(NR_contigs_Final$Final_Superkingdom == 'Bacteria', 'Bacteria', ifelse(NR_contigs_Final$Final_Superkingdom == 'Unclassified', 'Unclassified', ifelse(NR_contigs_Final$Final_Superkingdom == 'Archaea', 'Archaea', ifelse(NR_contigs_Final$Final_Class == "Caudoviricetes", "Prokaryotic Viruses", ifelse(NR_contigs_Final$Final_Family %in% euk_family, 'Eukaryotic Viruses', ifelse(NR_contigs_Final$Final_Family == "Unclassified", 'Unclassified Viruses', 'Prokaryotic Viruses'))))))))

### DATA-METADATA ASSOCIATION##
#Import metadata table (presented in supplementary table 1)
metadata <- read_delim("metadata_table.txt", "\t", escape_double = FALSE, col_names = TRUE, col_types = cols())

# Create summary table that will be progressively filled with the processed data
SummaryTable <- data.frame()

# Put all Viper magnitudes files in a separate folder. Then create a list with the file names.
list <- list.files()

# Import all magnitudes files from list, process the data file-by-file, write output in summary tables
#Import file and name the columns
for (i in list) {mag <- read_delim(i, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, col_types = cols())
colnames(mag) <- c('qseqid', 'mapped_reads')

#Add taxonomy to the magnitudes table by merging OverviewSub with mag by qseqid
All <- merge(mag, NR_contigs_Final, by = "qseqid", all.x = TRUE, all.y = TRUE)
#Remove all lines without reads
All <- filter(All, !(mapped_reads == 0))

#Add total reads column and calculate contig abundance percentages.Then write sample total contig table to file.
All <- mutate(All, 'total_sample_reads' = sum(All$mapped_reads))
All <- mutate (All, 'Contig_abundance_perc' = (mapped_reads / total_sample_reads) * 100)
write.xlsx(All, paste0(i,"_allContigs_overview.xlsx"))

#Add sample name column (for merging all samples data)
All <- mutate(All, Sample = word(i, 1, sep = fixed(".")))

#Fill in taxonomy summary table
SummaryTable <- plyr::rbind.fill(SummaryTable, All)}

#Export taxonomy summary table
write.xlsx(SummaryTable, "Taxonomy_Summary_Table.xlsx")

#Make MasterTable by including metadata
MasterTable <- merge(SummaryTable, metadata, by = "Sample", all.x = TRUE, all.y = TRUE)


###COVERAGE & DEPTH###
##Before this step, prepare the horizontal coverage table by:
#a) Indexing the bam files from Viper Classify and recording stats for each sample, run:
#samtools index -@ threads sample.Classify.bam
#samtools coverage $line.Classify.bam | pigz -c -p 12 > $line.coverage.gz
#b)Creating the horizontal coverage overview table, run:
#zcat $line.coverage.gz | cut -f1,6 | sed 1d | sort -k1,1 | sed 's/\t/,/g' | \
#sed "1s/^/contig,$line\n/" | join --header -t "," stats.horizontal_coverage.csv - > temp.stats
#mv temp.stats stats.horizontal_coverage.csv
##Load and prepare the tables for integration
#Load the coverage table
Coverage <- read_delim("stats.horizontal_coverage.csv", ",", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE, col_types = cols())
#Rename the contig column to qseqid to facilitate merging with the mastertable. R
names(Coverage)[names(Coverage) == 'contig'] <- 'qseqid'
#Transform to coverage table to long format
Cov_long <- melt(setDT(Coverage), id.vars = c("qseqid"), variable.name = "sample")
#Remove all the 0 coverage rows
Cov_long <- filter(Cov_long, !Cov_long$value == "0")
#Rename the value column to coverage
names(Cov_long)[names(Cov_long) == 'value'] <- 'coverage(perc)'

##Merge Coverage table with Master table by sample keeping the qseqids of the Master table
#Make sample vector to go sample-by-sample
sam <- unique(MasterTable$Sample)
#Create an empty data frame that will be filled sample-by-sample with the merged Master and Coverage tables
MasterTable_final <- data.frame()

for (p in sam) {temp1 <- filter(MasterTable, MasterTable$Sample == p)
temp2 <- filter(Cov_long, Cov_long$sample == p)
temp3 <- merge(temp1, temp2, by = "qseqid", all = TRUE)
MasterTable_final <- plyr::rbind.fill(MasterTable_final, temp3)}

#Remove the additional "samples" column from the table
MasterTable_final <- subset(MasterTable_final, select = -c(sample))


###FIGURES###
##Heatmap: Most abundant viral species by sample (log of mapped reads) - Coverage 70% - Relative abundance 0.005%
#Create a subset of the master table that includes only the eukaryotic viruses
euk_vir_master <- filter(MasterTable_final, Findings == "Eukaryotic Viruses")

#Filter the eukaryotic virus subset to include only results with coverage >= 70%
euk_vir_master70 <- subset(euk_vir_master, euk_vir_master$`coverage(perc)`>=70)
#Create a vector with all the unique eukaryotic viral species
vir <- unique(euk_vir_master70$Final_Species) 

HeatMap70 <- data.frame() #Create an empty data frame that will be filled with the abundance data for each eukaryotic virus for each sample

for (j in vir) {temp1 <- filter(euk_vir_master70, euk_vir_master70$Final_Species == j) #Create a temporary data frame for each unique virus
sam <- unique(temp1$Sample) #Create a vector for the samples that had reads for the unique virus
for (q in sam) {temp2 <- filter(temp1, temp1$Sample == q) #Create a temporary data frame of each sample for each unique virus
temp2 <- mutate(temp2, 'total_eukvir_reads' = sum(temp2$mapped_reads)) #Add all the reads mapped to each unique virus, regardless of the NR contig, and add them in a separate column
temp2 <- mutate (temp2, 'eukvir_abundance_perc' = (total_eukvir_reads / total_sample_reads) * 100) #Calculate the abundance of each virus and add it in a separate column
HeatMap70 <- plyr::rbind.fill(HeatMap70, temp2[1,])}} #For each sample for each virus write only the first row of the temporary data frame, the necessary columns are the same for all the row

#Filter out the reads that could not be classified at the species level
HeatMap70 <- filter(HeatMap70, !HeatMap70$Final_Species == "Unclassified")
#Calculate the log of viral reads for better depiction
HeatMap70 <- mutate(HeatMap70, 'log_eukvir_reads' = log2(HeatMap70$total_eukvir_reads))
#Set abundance cut-off for figure (0.005 used)
int <- subset(HeatMap70, HeatMap70$eukvir_abundance_perc>=0.005) 

#Fill in the low abundance cases for the the economically important viruses
#Create low abundance dataset 
out <- subset(HeatMap70, HeatMap70$eukvir_abundance_perc<0.005)
#Create low abundance dataframes for the economically important viruses (no additions needed for ToTV, CMV)
OutPepMV <- subset(out, out$Final_Species == "Pepino mosaic virus")
OutToBRFV <- subset(out, out$Final_Species == "Tomato brown rugose fruit virus")
OutSTV <- subset(out, out$Final_Species == "Amalgavirus lycopersici")
OutToFBV <- subset(out, out$Final_Species == "Blunervirus solani")
#Removed S47 from ToBRFV because it doesn't have >70% reference genome coverage
OutToBRFV <- subset(OutToBRFV, !OutToBRFV$SampleCode == "S47")
#Bind the viral low abundance dataframes to the heatmap dataframe
int <- rbind(int, OutPepMV, OutToBRFV, OutSTV, OutToFBV)
#Transform the data to the long data frame necessary for the heatmap
long <- subset(int, select = c("Final_Species", "SampleCode", "log_eukvir_reads"))

#Make the plot
#Sample type first then chronological order
a <- factor(long$SampleCode, levels = c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S11", "S12", "S18", "S19", "S20", "S21", "S23", "S25", "S26", "S27", "S28", "S29", "S30", "S31", "S32", "S33", "S34", "S35", "S36", "S37", "S38", "S39", "S40", "S41", "S42", "S43", "S44", "S45", "S46", "S47", "S48", "S49", "S50", "S51", "S52", "S53", "S54", "S55", "S56", "S57", "S58", "S59", "S60", "S62", "S63", "S64", "S65", "S67", "S68", "S69", "S70", "S71", "S72", "S73", "S74", "S75", "S76", "S77", "S79", "S80", "S81", "S82", "S83", "S84", "S85", "S86", "S87", "S88", "S91", "S92", "S93", "S94", "S95", "S97", "S99", "S9", "S10", "S13", "S14", "S22", "S66", "S78", "S89", "S90", "S96", "S98", "S100", "S101", "S15", "S16", "S17", "S24", "S61"))
#Host order
b <- factor(long$Final_Species, levels = c("Pepino mosaic virus", "Tomato brown rugose fruit virus", "Amalgavirus lycopersici", "Blunervirus solani", "Torradovirus lycopersici", "Cucumber mosaic virus", "Tomato mild mottle virus", "Potyvirus trompetae", "Botrytis cinerea umbra-like virus 1", "Sclerotinia sclerotiorum umbra-like virus 3-WX2", "Penicillimonavirus betaplasmoparae", "Penicillimonavirus etaplasmoparae", "Penicillimonavirus betapenicillii", "Penicillimonavirus kilnbarnense", "Sclerotinia sclerotiorum partitivirus 2", "Botryotinia fuckeliana partitivirus 1", "Sclerotinia sclerotiorum partitivirus 3", "Penicillium stoloniferum virus F", "Alphachrysovirus penicillii", "Erysiphe necator associated chrysovirus 2", "Erysiphe necator associated fusarivirus 3", "Sclerotinia sclerotiorum hypovirus 7", "Thetahypovirus sclerotiniae", "Botrytis virus F", "Monilinia barnavirus A", "Penicillum brevicompactum polymycovirus 1", "Cladosporium ramotenellum polymycovirus 1", "Xian Totiv tick virus 2", "Luoyang Totiv tick virus 2", "Totiviridae sp.", "Geotrichum candidum totivirus 1", "Trichoderma harzianum orthocurvulavirus 1", "Gemycircularvirus derva1", "Genomoviridae sp.", "Alphamesonivirus cavallyense", "Botrytis cinerea alpha-like virus 1", "Guyuan tick virus 1", "Psittaciform parvoviridae sp.", "Sichuan mosquito circovirus 3", "Tick-associated circular virus-6", "Fringilla montifringilla Circoviridae sp.", "Circoviridae sp.", "Solendovirus venanicotianae"))

ggplot(long, aes(x=a, y=b, fill=log_eukvir_reads))+
  geom_tile()+
  scale_x_discrete("Sample")+
  scale_y_discrete("Species", limits=rev, labels= c("Tobacco vein clearing virus", "Circoviridae sp.", "Fringilla montifringilla Circoviridae sp.", "Tick-associated circular virus-6", "Sichuan mosquito circovirus 3", "Psittaciform parvoviridae sp.", "Guyuan tick virus 1", "Botrytis cinerea alpha-like virus 1", "Alphamesonivirus 1 Ngewotan", "Genomoviridae sp.", "Tick-associated genomovirus 1", "Trichoderma harzianum orthocurvulavirus 1", "Geotrichum candidum totivirus 1", "Totiviridae sp.", "Luoyang Totiv tick virus 2", "Xian Totiv tick virus 2", "Cladosporium ramotenellum polymycovirus 1", "Penicillum brevicompactum polymycovirus 1", "Monilinia barnavirus A", "Botrytis virus F", "Sclerotinia sclerotiorum hypovirus 2", "Sclerotinia sclerotiorum hypovirus 7", "Erysiphe necator associated fusarivirus 3", "Erysiphe necator associated chrysovirus 2", "Penicillium raistrickii chrysovirus 1", "Penicillium stoloniferum virus F", "Sclerotinia sclerotiorum partitivirus 3", "Botryotinia fuckeliana partitivirus 1", "Sclerotinia sclerotiorum partitivirus 2", "Kiln Barn virus", "Penicillium glabrum negative-stranded RNA virus 1", "Plasmopara viticola lesion associated mononegaambi virus 9", "Plasmopara viticola lesion associated mononegaambi virus 2", "Sclerotinia sclerotiorum umbra-like virus 3-WX2", "Botrytis cinerea umbra-like virus 1", "Colombian datura virus", "Tomato mild mottle virus", "Cucumber mosaic virus", "Tomato torrado virus", "Tomato fruit blotch virus", "Southern tomato virus", "Tomato brown rugose fruit virus", "Pepino mosaic virus"))+
  scale_fill_distiller("Reads(log2)", palette ="Greens", direction = 1)+
  theme_bw()+ 
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())+
  theme(axis.text.y = element_text(face='italic', size = 14), axis.title.y = element_blank())+
  theme(legend.position="none")
#In Photoshop add: i)icons for host, ii)family information

##PepMV genotypes supplementary heatmap
#Filter the eukaryotic virus subset to include only results with coverage >= 70%
pepmv_70 <- subset(euk_vir_master70, euk_vir_master70$Final_Species == "Pepino mosaic virus")

CH2_70 <- subset(pepmv_70, pepmv_70$sseqid == "AZQ25016.1")
LP_70 <- subset(pepmv_70, pepmv_70$sseqid == "AVG19361.1")
EU_70 <- subset(pepmv_70, pepmv_70$sseqid == "WBG54278.1")
#Identify common contigs in EU and LP lists. Blast to clarify if finding is EU or LP. Check reference genome coverage and if contigs need to be removed.
#Make sample vectors with checked samples
true_ch2 <- c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19", "S20", "S21", "S22", "S23", "S24", "S25", "S26", "S27", "S28", "S29", "S30", "S31", "S32", "S33", "S34", "S35", "S36", "S37", "S38", "S39", "S40", "S41", "S42", "S43", "S44", "S45", "S46", "S47", "S48", "S49", "S50", "S51", "S52", "S53", "S54", "S55", "S56", "S57", "S58", "S59", "S60", "S61", "S62", "S63", "S64", "S65", "S66", "S67", "S68", "S69", "S70", "S71", "S72", "S73", "S74", "S75", "S76", "S77", "S78", "S79", "S80", "S81", "S82", "S83", "S84", "S85", "S86", "S87", "S88", "S89", "S90", "S91", "S92", "S93", "S94", "S95", "S96", "S97", "S98", "S99", "S100", "S101")
true_lp <- c("S34", "S35", "S36", "S53", "S54", "S57", "S63", "S69", "S70", "S78", "S80", "S87", "S88", "S90", "S93", "S95", "S96", "S97")
true_eu <- c("S30", "S31", "S41", "S48", "S60", "S61", "S68", "S74", "S8", "S86", "S91", "S92")

CH2_70 <- subset(CH2_70, CH2_70$SampleCode %in% true_ch2)
LP_70 <- subset(LP_70, LP_70$SampleCode %in% true_lp)
EU_70 <- subset(EU_70, EU_70$SampleCode %in% true_eu)

int_sup <- rbind(CH2_70, LP_70, EU_70)
long_sup <- subset(int_sup, select = c("sseqid", "SampleCode", "log_mapped_reads"))
#Sample type first then chronological order
a <- factor(long_sup$SampleCode, levels = c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S11", "S12", "S18", "S19", "S20", "S21", "S23", "S25", "S26", "S27", "S28", "S29", "S30", "S31", "S32", "S33", "S34", "S35", "S36", "S37", "S38", "S39", "S40", "S41", "S42", "S43", "S44", "S45", "S46", "S47", "S48", "S49", "S50", "S51", "S52", "S53", "S54", "S55", "S56", "S57", "S58", "S59", "S60", "S62", "S63", "S64", "S65", "S67", "S68", "S69", "S70", "S71", "S72", "S73", "S74", "S75", "S76", "S77", "S79", "S80", "S81", "S82", "S83", "S84", "S85", "S86", "S87", "S88", "S91", "S92", "S93", "S94", "S95", "S97", "S99", "S9", "S10", "S13", "S14", "S22", "S66", "S78", "S89", "S90", "S96", "S98", "S100", "S101", "S15", "S16", "S17", "S24", "S61"))
#Genotype order
b <- factor(long_sup$sseqid, levels = c("AZQ25016.1", "WBG54278.1", "AVG19361.1"))

ggplot(long_sup, aes(x=a, y=b, fill=log_mapped_reads))+
  geom_tile()+
  scale_x_discrete("Sample")+
  scale_y_discrete("PepMV Genotype", limits=rev, labels= c("PepMV-LP", "PepMV-EU", "PepMV-CH2"))+
  scale_fill_distiller("Reads(log2)", palette ="Greens", direction = 1)+
  theme_bw()+ 
  theme(axis.text.x = element_text(angle = 90, size = 7), axis.title.x = element_blank())

##Findings overview for viruses of interest
#For viruses of interest cases below the abundance cut-offs checked via qPCR to include all interesting cases
Findings_overview <- read.xlsx("Findings_overview.xlsx")

concs2 <- c("Diagnostics - HTS" = "#CA6F1E", "HTS - Molecular confirmation" = "#2874A6", "HTS - No molecular confirmation" = "#AED6F1")

ggplot(Findings_overview)+
  geom_bar(aes(x = factor(Virus, levels = c('PepMV-CH2', 'PepMV-EU', 'PepMV-LP', 'ToBRFV', 'STV', 'ToTV', 'ToFBV', 'CMV')), y = number, fill = factor(Cases, levels = c('HTS - No molecular confirmation', 'HTS - Molecular confirmation', 'Diagnostics - HTS'))), position = "stack", stat = "identity")+
  geom_text(data=subset(Findings_overview, !Findings_overview$number == "0"), aes(x = factor(Virus, levels = c('PepMV-CH2', 'PepMV-EU', 'PepMV-LP', 'ToBRFV', 'STV', 'ToTV', 'ToFBV', 'CMV')), y = text, label = number), vjust = 1.5, colour = "black")+
  scale_y_continuous("Samples", limits = c(0, 102))+
  scale_x_discrete("Viral species of interest")+
  scale_fill_manual("Viral identification", values = concs2, limits = c('HTS - No molecular confirmation', 'HTS - Molecular confirmation', 'Diagnostics - HTS'))+
  theme_bw()+
  theme(legend.position = "bottom")+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14), legend.text = element_text(size = 14), legend.title = element_text(size = 14))

##Amount per year bar chart
year <- read.xlsx("year_overview.xlsx")

concs<- c("2017" = "#B9770E", "2018" = "#B9770E", "2019" = "#B9770E", "2020" = "#B9770E", "2021" = "#B9770E", "2022" = "#B9770E", "2023" = "#B9770E", "2024" = "#B9770E")

ggplot(year)+
  geom_bar(aes(x=factor(Year, levels = c("2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024")), y=samples, fill=factor(Year, levels = c("2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024"))), position = "dodge", stat = "identity")+
  geom_text(aes(x=factor(Year, levels = c("2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024")), y=text,  label = samples), vjust = 1.5, colour = "black")+
  scale_fill_manual(values = concs, limits = c("2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024"))+
  scale_y_continuous("Samples")+
  scale_x_discrete("Sampling year", limits = c("2017", "2018", "2019", "2020", "2021", "2022", "2023", "2024"))+
  theme_bw()+
  theme(legend.position = "none")

##Type of sample bar chart
Type <- read.xlsx("sample_type.xlsx")

concs<- c("Leaf" = "#7FB3D5", "Fruit" = "#7FB3D5", "Sepal" = "#7FB3D5")

ggplot(Type)+
  geom_bar(aes(x=Type, y=Samples, fill=Type), position = "dodge", stat = "identity")+
  geom_text(aes(x=Type, y=text, label = Samples), vjust = 1.5, colour = "black")+
  scale_fill_manual(values = concs, limits = c("Leaf", "Fruit", "Sepal"))+
  scale_y_continuous("Samples")+
  scale_x_discrete("Type of sample", limits = c("Leaf", "Fruit", "Sepal"))+
  theme_bw()+
  theme(legend.position = "none")

## STV association with symptoms
#Check STV positive findings and identify sample pairs that have at least one sample positive for STV
stv_70 <- subset(euk_vir_master70, euk_vir_master70$Final_Species == "Amalgavirus lycopersici")
stv_full <- subset(euk_vir_master, euk_vir_master$Final_Species == "Amalgavirus lycopersici")

##Make a table with all the complete pairs
pairs <- c("P1", "P3", "P4", "P5", "P6", "P7", "P10", "P11", "P19", "P21", "P22")

symp_cor_pairs <- data.frame()
for (k in pairs) {temp3 <- subset(stv_full, stv_full$Pair == k) #Create a temporary data frame that includes both samples of the pair
symp_cor_pairs <- plyr::rbind.fill(symp_cor_pairs, temp3)}

#For pair 22 the symptomatic sample does not have any STV reads so it doesn't appear in the table. Make a new row for the sample so that the figure is complete.
#Duplicate the row of the non-symptomatic sample of pair 22.
symp_cor_pairs <- rbind(symp_cor_pairs, symp_cor_pairs[rep(21,1), ])
#Change the appropriate values in the new row so that it represents the symptomatic sample of pair 22.
symp_cor_pairs[22, 2] = "S89" #Correct original sample name
symp_cor_pairs[22, 3] = 0 #Correct mapped reads
symp_cor_pairs[22, 37] = 0 #Correct contig abundance
symp_cor_pairs[22, 39] = "Yes" #Correct symptoms status
symp_cor_pairs[22, 46] = 0 #Correct coverage %

#Add column with log mapped reads
symp_cor_pairs <- mutate(symp_cor_pairs, 'log_mapped_reads' = log2(symp_cor_pairs$mapped_reads))
symp_cor_pairs[22, 47] = 0
write.xlsx(symp_cor_pairs, "symp_cor_pairs.xlsx")
#symp_cor_single <- mutate(symp_cor_single, 'log_mapped_reads' = log2(symp_cor_single$mapped_reads))
#write.xlsx(symp_cor_single, "symp_cor_single.xlsx")

#Make the plot
concs_pairs<- c("Yes" = "#DC7633", "No" = "#52BE80")
pairs <- factor(symp_cor_pairs$Pair, levels = c("P1", "P3", "P4", "P5", "P6", "P7", "P10", "P11", "P19", "P21", "P22"))

ggplot(symp_cor_pairs)+
  geom_bar(aes(x=pairs, y=log_mapped_reads, fill=factor(Symptoms, levels = c('Yes', 'No'))), position = "dodge", stat = "identity")+
  facet_grid(~factor(symp_cor_pairs$infection_number, levels = c("STV + PepMV-CH2", "STV + PepMV-CH2 + VoI")), scales = "free_x", space = "free")+
  scale_fill_manual(values = concs_pairs, limits = c("Yes", "No"))+
  scale_y_continuous("STV mapped reads (log2)")+
  scale_x_discrete("Sample Pair")+
  theme_bw()+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))+
  theme(strip.text.x = element_text(size = 12))+
  #guides(fill = guide_legend(title = "Symptoms"))
  guides(fill = "none")

ggplot()+
  geom_boxplot(symp_cor_pairs, mapping = aes(x=factor(Symptoms, levels = c('Yes', 'No')), y=log_mapped_reads, fill=factor(Symptoms, levels = c('Yes', 'No')))) +
  facet_grid(~factor(symp_cor_pairs$infection_number, levels = c("STV + PepMV-CH2", "STV + PepMV-CH2 + VoI")), scales = "free_x", space = "free")+
  scale_fill_manual(values = concs_pairs, limits = c("Yes", "No"))+
  theme_bw()+ 
  scale_x_discrete("Symptoms")+
  scale_y_continuous("Total STV mapped reads(log2)")+
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))+
  theme(strip.text.x = element_text(size = 12))+
  guides(fill = guide_legend(title = "Symptoms"))


