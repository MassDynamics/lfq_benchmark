import pandas as pd 
import numpy as np 
import os

import plotly.graph_objects as go
import copy
from itertools import combinations

import numpy as np
from IPython.display import display, HTML
import re 

from ..protein_table_loader import ProteinTableLoader
from ..plotter import Plotter
from ..confusion_matrix_calculator import ConfusionMatrixCalculator   

class IPRG2015Benchmarker():

    def __init__(self, experiment_home, mode):
        self.spike_amounts = pd.read_csv("../resources/IPRGSpikes.csv",
                            index_col = 0)
        self.spike_mapping = {"Ovalbumin": "P44015", 
                "Myoglobin": "P55752", 
                "Phosphorylase_b": "P44374", 
                "BetaGalactosidase": "P44983" , 
                "BSA": "P44683", 
                "CarbonicAnhydrase": "P55249"}
        
        self.protein_table = ProteinTableLoader().run(experiment_home,mode)
        print("total proteins in protein table: {}".format(self.protein_table.shape[0]))
        #for i in self.spike_mapping.values():
        #    assert i in self.protein_table.index, "protein {} not in protein_table".format(i)

        self.protein_table = self.replace_hidden_labels(self.protein_table, self.spike_mapping)
        
        #for i in self.spike_mapping.keys():
        #    assert i in self.protein_table.index, "protein {} not in protein_table".format(i)

        self.experiment_home = experiment_home
        self.mode = mode
    
    def run(self):
        n_proteins = self.get_n_proteins(self.protein_table)
        print("Number of Proteins Identified: ", n_proteins)

        print((self.protein_table["Q-value-ident"] > 0.01 ).sum())
        estLogFCs = self.get_estimated_logFCs() 
        trueLogFCs = self.get_theoretical_logFCs()
        q_val_table = self.get_q_val_table()

        first_comparison = self.protein_table.filter(regex = "logFC", axis =1).columns[0][5:]
        #fig = Plotter().volcano_plot(self.protein_table, first_comparison, True)
        fig = Plotter().scatter_log_fc_accuracy(estLogFCs, trueLogFCs, q_val_table)
        fig.write_html(os.path.join(self.experiment_home, "logFCdifferences.html"))

        print(ConfusionMatrixCalculator().run(self))

        return 
        
    def get_estimated_logFCs(self):
        protein_table = self.protein_table
        protein_table = protein_table[protein_table.ProteinId.apply(lambda x: x in self.spike_mapping.keys())]
        protein_table = protein_table.filter(regex = "(logFC)", axis = 1)
        return protein_table

    def get_theoretical_logFCs(self):
        spike_amounts = self.spike_amounts
        logFC_df = pd.DataFrame(index = spike_amounts.index)
        for comparison in list(combinations(spike_amounts.columns,2)):
            c1 = comparison[0]
            c2 = comparison[1]
            predlogFC = np.log2(spike_amounts.iloc[:,int(c1)-1]/spike_amounts.iloc[:,int(c2)-1])
            try:
                logFC_df["logFC sample"+c1+" - sample"+c2] = predlogFC
            except:
                logFC_df["logFC Sample"+c1+" - Sample"+c2] = predlogFC
        logFC_df.index.name = "ProteinId"
        return logFC_df
        
    def get_n_proteins(self, protein_table):
        if "Q-value-ident" in protein_table.columns:
            n_proteins = sum(protein_table["Q-value-ident"]<0.01)
            return n_proteins
        else:
            print("No identity q-val column, returning number of rows")
            return protein_table.shape[0]

    def get_n_true_positives(self, protein_table):
        protein_table = copy.deepcopy(protein_table)
        protein_table.index = protein_table.ProteinId
        protein_table = protein_table.filter(regex = "adj.P.Val", axis = 1)
        protein_table = protein_table[protein_table.apply(lambda x: x<0.05)]
        dif_exp_proteins = protein_table.fillna(0)
        true_positives = [protein for protein in dif_exp_proteins 
                            if protein in self.spike_mapping.values()]
        return len(false_positives)

    def get_q_val_table(self):
        protein_table = self.protein_table
        spike_mapping = self.spike_mapping
        protein_table = self.filter_for_iPRG_proteins(protein_table.copy(), spike_mapping)
        q_val_table = protein_table.filter(regex = "(adj.P.Val )|(ProteinId)", axis = 1)
        q_val_table = q_val_table.filter(regex = "adj.P.Val", axis = 1)

        return q_val_table

    def filter_for_iPRG_proteins(self, protein_table, spike_mapping):
        protein_table = protein_table[protein_table.ProteinId.apply(lambda x: x in spike_mapping.keys())]
        return protein_table

    def replace_hidden_labels(self, protein_table, spike_mapping):

        for key, val in spike_mapping.items():
            protein_table.index = protein_table.index.str.replace(val,key)

        if "ProteinId" in protein_table.columns:
            for key, val in spike_mapping.items():
                protein_table.ProteinId = protein_table.ProteinId.str.replace(val,key)
        
        return protein_table
