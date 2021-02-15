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


class HerStudyBenchmarker():

    def __init__(self, experiment_home, mode):

        self.experiment_home = experiment_home
        self.mode = mode
        self.protein_table = ProteinTableLoader().run(experiment_home,mode)
        self.protein_table = self.protein_table[~self.protein_table.isna().any(axis = 1)]

        if mode == "Perseus":
            self.protein_table.columns = ["Q-value-ident","adj.P.Val","logFC","ProteinId"]
        elif mode == "BYO":
            self.protein_table.columns = ["Q-value-ident","logFC","adj.P.Val","ProteinId"]

    def run():

        n_proteins = self.get_n_proteins(self.protein_table)
        print("Number of Proteins Identified: ", n_proteins)

        return 

    def get_n_proteins(self, protein_table):
        n_proteins = sum(protein_table["Q-value-ident"]<0.01)
        return n_proteins

    def get_estimated_logFCs(self):
        cols_regex = "(logFC)"
        protein_table = self.protein_table.filter(regex = cols_regex , axis = 1)
        return protein_table

    def get_q_val_table(self):
        cols_regex = "(adj.P.Val)"
        protein_table = self.protein_table.filter(regex = cols_regex , axis = 1)
        return protein_table