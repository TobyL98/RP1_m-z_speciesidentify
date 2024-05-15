########################
# rules_cigar_scoresNQ

# this code determines which sites are common
# near the asparagines/ glutamines (N/Q) 
# that are deaminated
# looks 1,2 and 3 peptides before
# and 1, 2 and 3 petides after
# the N/Q that is deaminated
# Uses LC-MS/MS data
# These rules can then potentially be used to predict where
# N/Q are deaminated in theoretical in silico generated peptides
# Result is an excel file with the most common residues at each position


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

## Function - adds the before and after so we can complete the rules
def add_bef_aft(seq_bef, seq, seq_aft):
    return(seq_bef + seq + seq_aft)


##Finding patterns in the peptide sequence
#function will loops through all rows (i.e., each LC-MS/MS peptide fragment)
#ifoxidation count for site (P or K, N,Q is greater than 0)
#will find the aa target number away from site (P or K or N/Q)
#add them to a list
#takes in the LC-MS/MS dataframe, target_num(position that we are looking)
#num_of_targets - are there 1 or 2 targets
#target_num2 - second positions we are looking for if 2 targets
def before_O_site(df, site, target_num, num_of_targets, target_num2 = 0):
    #turns site into correct number for lookup
    if site == "P":
        site_num = 4
    elif site == "K":
        site_num = 2
    elif site == "N|Q":
        site_num = 1

    #loops through all rows in table
    target_list = []
    for index, row in df.iterrows():
        column = "deam_count"
        if row[column] > 0:
            #pulls out sequence and modification cigar strip
            seq = row["pep_seq"]
            mod = row['pep_var_mod_pos']
            mod = re.sub(r"^0.|.0$", "", mod) #cleans modiftcaion cigar strip
            index_list = []
            index = 0
            #goes through modification strips(e.g., 000400020)
            #if number matches number we want will put index in list
            for num in mod:
                if int(num) == site_num:
                    index_list.append(index)
                index += 1
            #uses list and target position to generate  target_list
            for index in index_list:
                if index +target_num < len(seq) and num_of_targets == 1:
                    target = seq[index + target_num]
                    target_list.append(target)
                #if two targets will create consensus
                elif index + target_num < len(seq) and num_of_targets == 2:
                    target = seq[index + target_num] + "X" + site + seq[index + target_num2]
                    target_list.append(target)


    
    #list of target site is used to create a dictionary with aa as the key
    #and the number of that aa in the list as the value (e.g., G: 100)
    target_dict = {}
    for x in target_list:
        if x in target_dict.keys():
            target_dict[x] = target_dict.get(x) + 1
        else:
            target_dict[x] = 0
            target_dict[x] = target_dict.get(x) + 1
    
    #converts dictionary to df
    pep_bef_df = pd.DataFrame.from_dict(target_dict, orient = 'index', columns = ["Count"])
    #calculates the percentages of each aa out of total
    total_count = pep_bef_df["Count"].sum()
    pep_bef_df["Percentage"] = round(pep_bef_df["Count"] / total_count * 100, 2)
    pep_bef_df = pep_bef_df.sort_values(by = ["Count"], ascending= False)
    return pep_bef_df


                        

    


#load in data
#load in all csv files in PTM_rules folder
csv_files = glob.glob('C:/Users/tobyl/OneDrive - The University of Manchester/Bioinformatics Masters/Research project 1/Git_repositories/RP1_m-z_speciesidentify/PTM_rules/LCMSMS/*.csv')

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


#creates dictionaries to store df's generated
pep_aft1_df_dict = {}
pep_aft2_df_dict = {}
pep_aft3_df_dict = {}
pep_bef1_df_dict = {}
pep_bef2_df_dict = {}
pep_bef3_df_dict = {}
pep_cons_df_dict = {}

#loops through the pep confidence scores
#and creates a df forc each confidence score
#assigns df to correct dictionary
pepScore_list = [0, 10, 20, 30, 40, 50]
for score in pepScore_list:
    #filters df by score
    df = df.loc[df["pep_score"] >= score]

    #use apply on function to count number of hydroxylations (P and K)
    #and number of demidations (N/Q)
    #from LC-MS/MS data
    df["Phyd_count"] = df["pep_var_mod"].apply(mod_count, res = "P")
    df["Khyd_count"] = df["pep_var_mod"].apply(mod_count, res = "K")
    df["deam_count"] = df["pep_var_mod"].apply(mod_count, res = "NQ")


    #turns before_O_site function
    #site is the where modification occurs
    #target_num is position away from site being targeted (+1 is one before site)
    pep_aft1_df_p = before_O_site(df, site = "N|Q", target_num = +1, num_of_targets = 1)
    score_heading = "pep_score = {}".format(score)
    #dictionary has the pep_score as a key and the df generated as a value
    pep_aft1_df_dict[score_heading] = pep_aft1_df_p

    pep_aft2_df_p = before_O_site(df, site = "N|Q", target_num = +2, num_of_targets = 1)
    score_heading = "pep_score = {}".format(score)
    pep_aft2_df_dict[score_heading] = pep_aft2_df_p

    pep_aft3_df_p = before_O_site(df, site = "N|Q", target_num = +3, num_of_targets = 1)
    score_heading = "pep_score = {}".format(score)
    pep_aft3_df_dict[score_heading] = pep_aft3_df_p

    pep_bef1_df_p = before_O_site(df, site = "N|Q", target_num = -1, num_of_targets = 1)
    score_heading = "pep_score = {}".format(score)
    pep_bef1_df_dict[score_heading] = pep_bef1_df_p

    pep_bef2_df_p = before_O_site(df, site = "N|Q", target_num = -2, num_of_targets = 1)
    score_heading = "pep_score = {}".format(score)
    pep_bef2_df_dict[score_heading] = pep_bef2_df_p

    pep_bef3_df_p = before_O_site(df, site = "N|Q", target_num = -3, num_of_targets = 1)
    score_heading = "pep_score = {}".format(score)
    pep_bef3_df_dict[score_heading] = pep_bef3_df_p
    

    #pep_cons_df_p = before_O_site(df, site = "N|Q", target_num = -2, target_num2 = 1, num_of_targets = 2)
    #score_heading = "pep_score = {}".format(score)
    #pep_cons_df_dict[score_heading] = pep_cons_df_p

#for each target position generates a single table with scores as the headings
pep_aft1_df_all = pd.concat(pep_aft1_df_dict.values(), keys= pep_aft1_df_dict.keys(), axis = 1)
pep_aft2_df_all = pd.concat(pep_aft2_df_dict.values(), keys= pep_aft2_df_dict.keys(), axis = 1)
pep_aft3_df_all = pd.concat(pep_aft3_df_dict.values(), keys= pep_aft3_df_dict.keys(), axis = 1)
pep_bef1_df_all = pd.concat(pep_bef1_df_dict.values(), keys= pep_bef1_df_dict.keys(), axis = 1)
pep_bef2_df_all = pd.concat(pep_bef2_df_dict.values(), keys= pep_bef2_df_dict.keys(), axis = 1)
pep_bef3_df_all = pd.concat(pep_bef3_df_dict.values(), keys= pep_bef3_df_dict.keys(), axis = 1)

#pep_cons_df_all = pd.concat(pep_cons_df_dict.values(), keys = pep_cons_df_dict.keys(), axis = 1)
#print(pep_cons_df_all)

#writes all the tables to an excel sheet
with pd.ExcelWriter("NQ_rules.xlsx") as writer:
    pep_aft1_df_all.to_excel(writer, sheet_name = "1AfterNQ")
    pep_aft2_df_all.to_excel(writer, sheet_name = "2AfterNQ")
    pep_aft3_df_all.to_excel(writer, sheet_name = "3AfterNQ")
    pep_bef1_df_all.to_excel(writer, sheet_name = "1BeforeNQ")
    pep_bef2_df_all.to_excel(writer, sheet_name = "2BeforeNQ")
    pep_bef3_df_all.to_excel(writer, sheet_name = "3BeforeNQ")