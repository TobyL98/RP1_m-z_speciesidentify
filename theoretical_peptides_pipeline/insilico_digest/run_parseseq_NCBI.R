#######################
## run_parseseq_NCBI.R
#######################

# take dataframe  from in_silico_res_NCBI
# produces peptides generated from that for each species
# using parse_seq_M function
# this can then be merged with LCMSMS data in integrate_NCBI.py
# Code adapted from bacollite, please go there for further information or help

library("dplyr")
library("stringr")
library("bacollite")


##adding the source for the functions used
# use full filepath for subprocess run
source("C:/Users/tobyl/OneDrive - The University of Manchester/Bioinformatics Masters/Research project 1/Git_repositories/RP1_m-z_speciesidentify/theoretical_peptides_pipeline/insilico_digest/parseseq_M.R")
source("C:/Users/tobyl/OneDrive - The University of Manchester/Bioinformatics Masters/Research project 1/Git_repositories/RP1_m-z_speciesidentify/theoretical_peptides_pipeline/insilico_digest/mass_iso_M.R")
source("C:/Users/tobyl/OneDrive - The University of Manchester/Bioinformatics Masters/Research project 1/Git_repositories/RP1_m-z_speciesidentify/theoretical_peptides_pipeline/insilico_digest/RcppExports.R")

###Read in csv file
seq_df <- read.csv("C:/Users/tobyl/OneDrive - The University of Manchester/Bioinformatics Masters/Research project 1/Git_repositories/RP1_m-z_speciesidentify/theoretical_peptides_pipeline/Sequences/sequences_taxon_NCBI.csv", sep = ",")
print("Creating theoretical peptides")
print("Slow part of code. Outputs are below:")

#############
##Function
#############

#loops through all the rows in dataframe
#uses parse.seq to calculate peptide fragments and masses for each fragment
#results is dataframe with masses of peptide fragments
#for every sequence in the original dataframe
mass_loop <- function(df){
  for (row in 1:nrow(df)){
    sequence <- df[row, "sequence"] #takes sequence
    organism <- df[row, 
                   c("GENUS", "SPECIES")] #take organism
    
    csv_name <- df[row, "SPECIES"]
    csv_name <- gsub(" ", "_", csv_name)
    csv_name <- paste("C:/Users/tobyl/OneDrive - The University of Manchester/Bioinformatics Masters/Research project 1/Git_repositories/RP1_m-z_speciesidentify/theoretical_peptides_pipeline/insilico_digest/in_silico_res_NCBI/"
                      , csv_name, ".csv", sep = "")
    print(csv_name)
    #calculates peptide fragments and masses
    pepmass_df <- parse.seq(sequence, max.missed.cleaves = 1)
    #adds organism name as identifier
    pepmass_df <- cbind(pepmass_df, organism)
    write.csv(pepmass_df, file = csv_name)
    
  }

}

mass_loop(seq_df)
print("######################################")



