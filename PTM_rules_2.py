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


## Counting total number of possible PTMs in sequence
def Tot_PTMs(seq, search):
    res_count = 0
    for res in seq:
        if re.match(search, res):
            res_count += 1
    return(res_count)        


##Loading modules
import pandas as pd
import matplotlib as plt
import re

#load in data
df = pd.read_csv("LCMSMS_data_cattle.csv", sep = ",")

#filter dataframne to keep three columns
df = df[["pep_seq", "pep_var_mod", "pep_var_mod_pos"]]

#use apply on function to count number of hydroxylations (P and K)
#from LC-MS/MS data
#and deamidations (NQ)
df["Phyd_count"] = df["pep_var_mod"].apply(mod_count, res = "P")
df["Khyd_count"] = df["pep_var_mod"].apply(mod_count, res = "K")
df["Tot_hyd_count"] = df["Khyd_count"] + df["Phyd_count"]
df["deam_count"] = df["pep_var_mod"].apply(mod_count, res = "NQ")

#use apply on function to calculate theoretically possible
#hydroxylations and deamidations from the sequence only
df["Max_hyd_count"] = df["pep_seq"].apply(Tot_PTMs, search = r"P|K")
df["Max_deam_count"] = df["pep_seq"].apply(Tot_PTMs, search = r"N|Q")


##Summarising the number of counts for each max count for hydroxylations
max = df["Max_hyd_count"].max() # max number of theoretically possible hydroxylations
app_perc_list = []

#for loop loops through all peptides have
#certain number of possible hydroxylations
for num in range(0, max):
    summary_df = df.loc[(df["Max_hyd_count"] == num)] #filters by number of possible hydroxylations
    #summarises counts of number of actual hydroxylations in Lc-MS/MS data
    app_perc = summary_df["Tot_hyd_count"].value_counts()
    app_perc = app_perc.rename(num)
    #appends the summarised count series to a list
    app_perc_list.append(app_perc)

# converts all series in the list into a df  
hyd_rules_df = pd.concat(app_perc_list, axis = 1)
# converts to csv
hyd_rules_df.to_csv("hyd_counts.csv")



##Summarising the number of counts for each max count for deamidations
max = df["Max_deam_count"].max()
app_perc_list = []

#for loop loops through all peptides have
#certain number of possible deamidations
for num in range(0, max):
    summary_df = df.loc[(df["Max_deam_count"] == num)] #filters by number of possible deamidations
    #summarises counts of number of actual deamidations in Lc-MS/MS data
    app_perc = summary_df["deam_count"].value_counts()
    app_perc = app_perc.rename(num)
    #appends the summarised count series to a list
    app_perc_list.append(app_perc)
    
# converts all series in the list into a df    
deam_rules_df = pd.concat(app_perc_list, axis = 1)
# converts to csv
deam_rules_df.to_csv("deam_counts.csv")

