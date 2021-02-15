import os
import re
import pandas as pd
import numpy as np

class ProteinTableLoader():

    def run(self, experiment_home, mode):

        if mode == "MD":
            protein_table =  self.source_md_table(experiment_home)
        
        elif mode == "BYO":
            protein_table = self.source_byo_table(experiment_home)
        
        elif mode == "Perseus":
            if "PerseusResults" in os.listdir(experiment_home):
                protein_table = self.source_perseus_tables(os.path.join(experiment_home,
                                                                    "PerseusResults"))
            else:
                protein_table = self.source_perseus_tables(experiment_home)
        
        elif mode == "PD":
            protein_table = self.source_pd_table(experiment_home)
        else:
            print("Modes supported: MD+ Discovery, MD+ Discovery BYO, Perseus, Proteome Discoverer")

        self.check_protein_table_integrity(protein_table)

        if "Q-value-ident" in protein_table.columns:
            return protein_table[protein_table["Q-value-ident"]<0.01]
            
        else:
            print("No identity q-val column, all rows")
            return protein_table

    def source_md_table(self, experiment_home):

        if "transform" in os.listdir(experiment_home):
            protein_table = pd.read_csv(os.path.join(experiment_home,"transform",
                            "Protein_table.txt"), sep = "\t")
        else:
            protein_table = pd.read_csv(os.path.join(experiment_home,
                            "Protein_table.txt"), sep = "\t")
                            
        protein_table.ProteinId = protein_table.ProteinId.str.extract(r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})")[0]
       
        protein_table.index = protein_table.ProteinId
        
        # rename q-value column
        protein_table = protein_table.rename(columns=lambda x: re.sub('Protein q-value','Q-value-ident',x))
        
        # rename
        protein_table = protein_table.rename(columns=lambda x: re.sub('Sample','sample',x))
        
        # retain only required tables
        cols_regex = "(Q-value-ident)|(Difference)|(adj.P.Val )|(ProteinId)|(logFC)"
        protein_table = protein_table.filter(regex = cols_regex , axis = 1)
        return protein_table

    def source_byo_table(self,folder):
        
        if "transform" in os.listdir(folder):
            protein_table = pd.read_csv(os.path.join(folder,"transform",
                                "proteinGroups_quant.txt"), sep = "\t")
        else:
            protein_table = pd.read_csv(os.path.join(folder,
                                "proteinGroups_quant.txt"), sep = "\t")         
        print("Raw data shape: {}".format(protein_table.shape[0]))
        
        ids = protein_table["fasta headers"].apply(lambda x: str(x).split(" ")[0])

        #mask = [bool(1-i) for i in protein_table["majority protein ids"].str.contains("CON").to_list()]
        #protein_table = protein_table[mask]
        protein_table["ProteinId"] = ids.str.extract(r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})")[0]
       
        protein_table.index = protein_table.ProteinId

        # ensure that ident q-value isn't confused with adjusted values
        protein_table = protein_table.rename(columns=lambda x: re.sub('q-value','Q-value-ident',x))
        
        # retain only required tables
        cols_regex = "(Q-value-ident)|(Difference)|(adj.P.Val )|(ProteinId)|(logFC)"
        protein_table = protein_table.filter(regex = cols_regex , axis = 1)

        # correct column names for expression
        protein_table = protein_table.rename(columns=lambda x: re.sub('Welch\'s T-test ','',x))
        protein_table = protein_table.rename(columns=lambda x: re.sub('Sample','sample',x))
        protein_table = protein_table.rename(columns=lambda x: re.sub('_',' - ',x))

        return protein_table

    def source_perseus_tables(self, folder):

        comparison_files = [i for i in os.listdir(folder) if ".txt" in i]
        protein_table = pd.read_csv(os.path.join(folder,
                                comparison_files[0]), sep = "\t")

        #print("BYO Paper methods state that students T-tests are used.")
        #print("Please ensure that the correct columns are present here:")
        #print(protein_table.columns)

        print("Raw data shape: {}".format(protein_table.shape[0]))
        #protein_table = self.label_columns(protein_table, comparison_files[0])
        
        protein_table = self.format_perseus_column_names(protein_table)
        protein_table = protein_table.rename(columns=lambda x: self.fix_sample_names(x))


        if len(comparison_files) >1:
            for comparison in comparison_files[1:]:
                comparison_table = pd.read_csv(os.path.join(folder, comparison), sep = "\t")
                comparison_table = self.format_perseus_column_names(comparison_table)
                comparison_table = comparison_table.rename(columns=lambda x: self.fix_sample_names(x))
                protein_table = pd.merge(protein_table, comparison_table.drop("Q-value",1), 
                                        how = "outer",
                                        on = "Majority protein IDs")

        # remove junk rows
        #protein_table = protein_table.iloc[2:,:]

        # ensure that ident q-value isn't confused with adjusted values
        protein_table = protein_table.rename(columns=lambda x: re.sub('Q-value','Q-value-ident',x))
        ids = protein_table['Majority protein IDs'].str.extract(r"([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})")[0]
       
        cols_regex = "(Q-value-ident)|(logFC)|(adj.P.Val)"
        protein_table = protein_table.filter(regex = cols_regex , axis = 1)
        protein_table = self.format_perseus_column_names(protein_table)
        protein_table = protein_table.astype(float)
        protein_table["ProteinId"] = ids
        protein_table.index = protein_table.ProteinId

        #protein_table = protein_table.apply(lambda x: self.correct_q_value_scale_direction(x) if "adj.P.Val" in x.name else x)

        return protein_table

    def source_pd_table(self, experiment_home):

        protein_file = [i for i in os.listdir(experiment_home) if "proteins" in i][0]
        protein_table = pd.read_excel(os.path.join(experiment_home, protein_file))

        #get id
        protein_table["ProteinId"] = protein_table.Accession
        protein_table.index = protein_table.ProteinId
        
        #identification q-value
        protein_table = protein_table.rename(columns=lambda x: re.sub('Exp. q-value: Combined','Q-value-ident',x))
        
        #Log FC column: (pre-empt reversal of comparison sign)
        protein_table.columns = protein_table.columns.str.replace(r'Abundance Ratio: \(([0-9]+)\) \/ \(([0-9]+)\)',
                                                                    lambda x: "logFC sample{} - sample{}".format(x.group(2), x.group(1)))
        
        # log and reverse comparison
        cols = [col for col in protein_table.columns if "logFC" in col]
        for col in cols:
            protein_table[col] = -1*np.log2(protein_table[col])

        
        protein_table.columns = protein_table.columns.str.replace(r'Abundance Ratio Adj. P-Value: \(([0-9]+)\) \/ \(([0-9]+)\)',
                                                                    lambda x: "adj.P.Val sample{} - sample{}".format(x.group(2), x.group(1)))
        

        #filter to needed columns
        cols_regex = "(Q-value-ident)|(adj.P.Val )|(ProteinId)|(logFC)"
        protein_table = protein_table.filter(regex = cols_regex , axis = 1)

        #order columns
        protein_table = protein_table.reindex(sorted(protein_table.columns), axis=1)
        
        return protein_table
   
    def label_columns(self, protein_table, name):
        print(protein_table.columns)
        #print(name)
        c1 = name.split("_")[0]
        c2 = name.split("_")[1].split(".")[0]
        pattern = " sample" + c1 + " - sample" +c2
        #pattern = name.split(".")[0]
        protein_table = protein_table.rename(columns=lambda x: re.sub('-Log\(P-value\)','adj.P.Val' + pattern,x))
        protein_table = protein_table.rename(columns=lambda x: re.sub('Difference','Difference'+pattern,x))
        #protein_table = protein_table.rename(columns=lambda x: re.sub('Difference','Difference',x))

        return protein_table

    def correct_q_value_scale_direction(self, x):

        # if scale is correct:
        if max(x) < 1: # not on log scale
            print("adjusted p-values not on log scale")
            return x
        else: # on log scale
            print("exponentiate values")
            return 10**(-1*x) # fix log

    def fix_sample_names(self, name):
        
        if ("logFC" in name) or ("adj.P.Val" in name):

            name_prefix = name.split("_")[0]
            name_suffix = name.split("_")[1]

            name_prefix = name_prefix.split(" ")[0] + " sample"+name_prefix.split(" ")[1]
            name_suffix = " - sample" + name_suffix
            name = name_prefix + name_suffix
            print(name)
            return name
        
        else:
            return name

    def format_perseus_column_names(self, protein_table):
        protein_table = protein_table.rename(columns=lambda x: re.sub('C: ','',x))
        protein_table = protein_table.rename(columns=lambda x: re.sub('T: ','',x))
        protein_table = protein_table.rename(columns=lambda x: re.sub('N: ','',x))
        protein_table = protein_table.rename(columns=lambda x: re.sub('Difference','logFC',x))
        protein_table = protein_table.rename(columns=lambda x: re.sub('q-value','adj.P.Val',x))
        protein_table = protein_table.rename(columns=lambda x: re.sub('Student\'s T-test ','',x))
        protein_table = protein_table.rename(columns=lambda x: re.sub('Welch\'s T-test ','',x))
        return protein_table

    def check_protein_table_integrity(self, protein_table):

        # assert proteinId is present
        assert "ProteinId" in protein_table.columns, "ProteinId not in columns."
        
        # assert each logFC has a p-value to match the comparison:
        for logFC in protein_table.columns:
            if "logFC" in logFC:
                matching_cols = 0
                for pval in protein_table:
                    if "adj.P.Val" in pval:
                        if logFC.split("logFC")[1] == pval.split("adj.P.Val")[1]:
                            matching_cols = matching_cols + 1
                if not matching_cols != 0:
                   print("N matching cols for {} is 0".format(logFC))
                if not matching_cols <2:
                    print("N matching cols for {} are {}".format(logFC,
                                                                    matching_cols))


        # assert adj.p.val on 0 - 1
        for pval in protein_table:
            if "adj.P.Val" in pval:
                if not protein_table[pval].max() <= 1: 
                    print("{} may not be on scale 0-1, values more than 1".format(pval))
                if not protein_table[pval].min() >= 0:
                    print("{} may not be on scale 0-1, values less than 0".format(pval))
                #assert protein_table[pval].mean() > 0.5, "{} may not .format(pval)
                
        print("Protein table has passed integrity test.")
        return 
