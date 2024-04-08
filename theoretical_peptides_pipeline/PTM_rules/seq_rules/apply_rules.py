#################################
#CODE apply_rules.py
#################################

# Code applies the PTM rules calculated in
# for proline and lysine hydroxylation
# sees if it accounts for all of the total
# oxidations that were collected from LCMSMS data
# calculates percentage of oxidations accounted
# for by rules
# looks at each peptide in LCMSMS data to do this



##Loading modules
import pandas as pd
import matplotlib.pyplot as plt
import re
import glob

######################
## FUNCTIONS
######################

#This function uses the Pep-var_mod column
#to work out the number of PTMS in targeted residue (P, K, N/Q)
#that there is in each peptide fragment sequenced by LC-MS/MS
def mod_count(mods, res):
    mods = str(mods)
    
    if re.search(r";", mods):
        #splits if it has ; 
        # then seperates target residue modifications from other mods
        # example input: Oxidation (K); 2 Oxidation (P) 
        all = mods.split(";") 
        for x in all:
            if re.search(res, x):
                mods = x.strip(" ")

    #if it contains target residue
    if re.search(res, mods):
        mod_count = 1 #assumes it will be 1
        mod_num = mods.split(" ")[0] #splits to get number

        #if the number matches a number between 2 and 10
        #assigns that number to count, overwrites 1
        for num in range(2, 11):
            if mod_num == str(num):
                mod_count = mod_num

    #if it doesn't contain target residue count is 0        
    else:
        mod_count = 0
    return(int(mod_count))


def pred_p(seq):
    found = re.findall(r"(?=G[A-Z]PG)", seq)
    count = len(found)
    return count

def pred_k(seq):
    found = re.findall(r"(?=G[A-Z]KG)", seq)
    count = len(found)
    return count

#load in data
#load in all csv files in PTM_rules folder
csv_files = glob.glob('C:/Users/tobyl/OneDrive - The University of Manchester/Bioinformatics Masters/Research project 1/Git_repositories/RP1_m-z_speciesidentify/theoretical_peptides_pipeline/PTM_rules/LCMSMS/*.csv')

#loops through all CSV files
#creates dataframe with columns required
#adds df to list
df_list = []
for csv in csv_files:
    new_df = pd.read_csv(csv, sep = ",")
    new_df = new_df[["pep_res_before","pep_seq", "pep_res_after", "pep_var_mod", "pep_var_mod_pos","pep_score"]]
    df_list.append(new_df)

#combines all dataframes in the list into single df   
df = pd.concat(df_list)
df = df.loc[df["pep_score"] >= 30]
#use apply on function to count number of hydroxylations (P and K)
#from LC-MS/MS data
df["Phyd_count"] = df["pep_var_mod"].apply(mod_count, res = "P")
df["Khyd_count"] = df["pep_var_mod"].apply(mod_count, res = "K")
df["pep_seq_bef_aft"] = df["pep_res_before"] + df["pep_seq"] + df["pep_res_after"]

df["Pred_Phyd_count"] = df["pep_seq_bef_aft"].apply(pred_p)
df["Pred_Khyd_count"] = df["pep_seq_bef_aft"].apply(pred_k)
#do they match
countP = 0
countK = 0
for index, row in df.iterrows():
    if row["Phyd_count"] == row["Pred_Phyd_count"]:
        countP += 1
    if row["Khyd_count"] == row["Pred_Khyd_count"]:
        countK += 1
percentageP = (countP / df.shape[0]) * 100
percentageK = (countK/ df.shape[0]) * 100
print("The rule works {0} of times for Proline".format(percentageP))
print("The rule works {0} of times for Lysine".format(percentageK))

print(df.tail())
df.to_csv("rules_test.csv")




