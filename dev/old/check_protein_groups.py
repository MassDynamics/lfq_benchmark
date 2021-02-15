# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %%
from IPython import get_ipython

# %% [markdown]
# # Check Protein Groups

# %%
get_ipython().system(' pwd')


# %%
import pandas as pd


# %%
# check data once you remove 
ups  = pd.read_csv("../protein_groups/proteinGroups_UPS.txt", sep = "\t")
cols = ['Reverse','Only identified by site','Potential contaminant']
for col in cols:
    print(col)
    print(ups.shape)
    print(len(ups["Majority protein IDs"][ups["Majority protein IDs"].str.contains("ups")].unique()))
    print("apply filter")
    ups = ups[ups[col] != "+"]
    print(ups.shape)
    print(len(ups["Majority protein IDs"][ups["Majority protein IDs"].str.contains("ups")].unique()))
    print()
    


# %%
ups_proteins = pd.read_csv("../resources/UPS_benchmark.tsv", sep = "\t").iloc[:,0].to_list()
print(len(ups_proteins))
ups_proteins[0:3]


# %%
ups  = pd.read_csv("../protein_groups/proteinGroups_UPS.txt", sep = "\t")
ups["Majority protein IDs"][ups["Majority protein IDs"].str.contains("ups") & (ups["Potential contaminant"] == "+")]

# %% [markdown]
# # UPS Conclusion: 
# Of 48 UPS proteins, 2 reverse sequences were generated, adding to 50 in total. After filtering by filter groups, we would lose 2 reverse sequences (good) but 2 potential contaminents (ovalbumin and gelsin) that we want to keep. To retain themm we will manually edit the file (reloading to be sure). 
# The number of rows in the resulting data should be 2242. 

# %%
# check data once you remove 
ups  = pd.read_csv("../protein_groups/proteinGroups_UPS_remove_contaminent_flag_on_spiked_proteins.txt", sep = "\t")
cols = ['Reverse','Only identified by site','Potential contaminant']
for col in cols:
    print(col)
    print(ups.shape)
    print(len(ups["Majority protein IDs"][ups["Majority protein IDs"].str.contains("ups")].unique()))
    print("apply filter")
    ups = ups[ups[col] != "+"]
    print(ups.shape)
    print(len(ups["Majority protein IDs"][ups["Majority protein IDs"].str.contains("ups")].unique()))
    print()
    

# %% [markdown]
# ## Perfect!

# %%



