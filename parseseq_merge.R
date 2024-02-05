library("bacollite")
library("dplyr")

setwd("C:/Users/tobyl/OneDrive - The University of Manchester/Bioinformatics Masters/Research project 1/R work")
###Read in csv file
seq_df <- read.csv("seq.csv", sep = ",")

#rename column to organism
seq_df <- seq_df %>%
  rename(
    Organism = X
  )

##Function
#loops through all the rows in dataframe
#uses parse.seq to calculate peptide fragments and masses for each fragment
#results is dataframe with masses of peptide fragments
#for every sequence in the original dataframe
mass_loop <- function(df){
  for (row in 1:nrow(df)){
    sequence <- df[row, "sequence"] #takes sequence
    organism <- df[row, "Organism"] #take organism
    
    if (row == 1){
      #calculates peptide fragments and masses
      pepmass_df <- parse.seq(sequence, max.missed.cleaves = 1)
      #adds organism name as identifier
      pepmass_df <- cbind(pepmass_df, organism)
    }
    else {
      #calculates peptide fragments and masses
      pepmass_df2 <- parse.seq(sequence, max.missed.cleaves = 1)
      #adds organism name as identifier
      pepmass_df2 <- cbind(pepmass_df2, organism)
      
      #merges this organism df with rest of df
      pepmass_df <- rbind(pepmass_df, pepmass_df2)
    }
  }
  return(pepmass_df)
}

pep_df <- mass_loop(seq_df)

