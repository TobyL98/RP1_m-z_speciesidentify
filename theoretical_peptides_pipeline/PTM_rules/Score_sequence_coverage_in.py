##### Calculating the sequence coverage for different LC-MS/MS scores

import pandas as pd
import matplotlib.pyplot as plt
import re

df = pd.read_csv("results.csv", sep = ',')

#####Taking out the protein sequence
#prot_sequence = df.loc[1, "prot_seq"]
total_coverage = 2057

score_list = [0, 10, 20, 30, 40, 50]
perc_coverage_list = []
#loops through the peptide_confidence score
#will then filter by confidence score
#and calculate the percentage coverage from that
for score in score_list:
    #calculating coverage at different peptide_scores
    df = df[["pep_start", "pep_end", "pep_score", "pep_seq"]]
    df = df.loc[df["pep_score"] >= score]

    #sorting be pep_score value
    df = df.sort_values("pep_score", ascending= False)
    #dropping duplicates
    #have to drop any duplicates with overlapping sequence coverage
    df = df.drop_duplicates(["pep_start", "pep_end", "pep_seq"])
    df = df.drop_duplicates(["pep_start"])
    df = df.drop_duplicates(["pep_end"]).sort_index()
    df = df.sort_values(["pep_start"])
    print(df.head())
    #calculate the length of the peptide
    df["pep_len"] = (df["pep_end"] - df["pep_start"]) + 1

    #calculating the percentage coverage of the LC-MS/MS peptides
    pep_coverage = df["pep_len"].sum()
    perc_coverage = (pep_coverage / total_coverage) * 100
    perc_coverage_list.append(perc_coverage)

 #plotting the results using matplotlib   
plt.plot(score_list, perc_coverage_list)
plt.xlabel("LC-MS/MS Peptide Score", fontsize = 15)
plt.ylabel("Sequence Coverage (%)", fontsize = 15)
plt.show()
#plt.savefig("score_sequenceCoverage.png")

