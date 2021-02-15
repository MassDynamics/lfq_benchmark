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


class UPSBenchmarker():

    def __init__(self, experiment_home):

        self.experiment_home = experiment_home
        self.spike_amounts = self.get_usp_data()
        self.mode=mode
        self.protein_table =  ProteinTableLoader().run(experiment_home,mode)

    def run(self):

        trueLogFCs = self.get_theoretical_logFCs()
        estLogFCs = self.get_estimated_logFCs()
        q_val_table = self.get_q_val_table()


        n_proteins = self.get_n_proteins(self.protein_table)
        print("Number of Proteins Identified: ", n_proteins)
        
        first_comparison = self.protein_table.filter(regex = "logFC", axis =1).columns[0][5:]
        #fig = Plotter().volcano_plot(self.protein_table, first_comparison, True)
        fig = Plotter().scatter_log_fc_accuracy(estLogFCs, trueLogFCs, q_val_table)
        fig.write_html(os.path.join(self.experiment_home, "logFCdifferences.html"))

        print(ConfusionMatrixCalculator().run(self))
    
        return 
        
    def get_usp_data(self):
        usp_data = pd.read_csv("../resources/UPS_benchmark.tsv", sep = "\t")
        usp_data = usp_data.iloc[:,:3]
        usp_data.columns = ["ProteinId","UPS1","UPS2"]
        usp_data.index = usp_data.ProteinId
        usp_data["UPS1"] = usp_data.UPS1.str.replace(",","").astype(float)
        usp_data["UPS2"] = usp_data.UPS2.str.replace(",","").astype(float)
        usp_data = usp_data.drop("ProteinId", axis = 1)

        return usp_data

    def get_theoretical_logFCs(self):
        protein_table = self.get_usp_data()
        protein_table["logFC UPS1 - UPS2"] = np.log2(protein_table.UPS1.div(protein_table.UPS2))
        protein_table = protein_table[["logFC UPS1 - UPS2"]]
        return protein_table

    def get_q_val_table(self):
        protein_table = self.protein_table
        usp_proteins = [i for i in self.get_usp_data().index.to_list() if i in protein_table.index]
        protein_table = protein_table.loc[usp_proteins]
        cols_regex = "(adj.P.Val )"
        protein_table = protein_table.filter(regex = cols_regex , axis = 1)
        return protein_table

    def get_naive_estimated_logFCs(self):
        protein_table = self.protein_table
        in_table = [i for i in self.get_usp_data().index.to_list() if i in protein_table.index]
        out_table = [i for i in self.get_usp_data().index.to_list() if i not in protein_table.index]

        if out_table:
            print('The following spiked proteins were not identified:')
            for i in out_table:
                print(i)
            print(len(out_table), " proteins were missing in total")
        protein_table = protein_table.loc[in_table]
        cols_regex = "(logFC)"
        protein_table = protein_table.filter(regex = cols_regex , axis = 1)
        return protein_table

    def get_estimated_logFCs(self):

        return self.get_naive_estimated_logFCs()

    def get_n_proteins(self, protein_table):
        if "Q-value-ident" in protein_table.columns:
            n_proteins = sum(protein_table["Q-value-ident"]<0.01)
            return n_proteins
        else:
            print("No identity q-val column, returning number of rows")
            return protein_table.shape[0]

    def get_adjusted_estimated_logFCs(self):

        logFCs = self.get_theoretical_logFCs()
        true_0_index = logFCs[logFCs.iloc[:,0] == 0].index
        false_0_mean = self.get_naive_estimated_logFCs().loc[true_0_index].mean().values[0]
        print("Dynamic Benchmark Dataset estimated false 0 mean: ", false_0_mean)
        adj_LogFCs = self.get_naive_estimated_logFCs().apply(lambda x: x - false_0_mean, axis = 1)

        return adj_LogFCs
