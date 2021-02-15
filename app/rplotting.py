from functools import partial
from rpy2 import robjects
from rpy2.ipython import html
import rpy2.ipython.html
from rpy2.robjects.packages import importr
import rpy2.robjects.packages as rpackages
import os
from IPython.display import Image
import numpy as np
from scipy.stats import pearsonr
import pandas as pd


# required for running R plotting stuff
base = importr('base')
utils = importr('utils')

def install_r_packages(packnames):
    from rpy2.robjects.vectors import StrVector
    names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
    utils.chooseCRANmirror(ind=1)
    if len(names_to_install) > 0:
        utils.install_packages(StrVector(names_to_install))

install_r_packages(["ggplot2","hexbin","RColorBrewer"])


ggplot2 = rpackages.importr('ggplot2')
hexbin = rpackages.importr('hexbin')
brewer = rpackages.importr('RColorBrewer')

class RPlotting():

    def hexbin(self, table, save_file = False, 
                xlab = "", ylab = "",
                lims = "c(-13,13)", 
                nbins = 30):

        # set default name
        if not save_file:
            name = "rplot.png"
        else:
            name = save_file

        # write table
        table.to_csv("tmp.csv")

        #print r2 
        score = pearsonr(table.x,table.y)[0]
        print('Predictions Hexbin: pearson corr = ' + str(round(score,3)))
        robjects.r(
        '''
        # get data
        table = read.csv("tmp.csv")

        #make plot
        bin = hexbin(table$x, table$y, xbnds = {}, ybnds = {}) 
        my_colors=colorRampPalette(rev(brewer.pal(11,'Spectral')))
        png("rplot.png") 
        plot(bin, main="" , colramp=my_colors, 
                xlab = {}, 
                ylab = {}) 
        dev.off()
        '''.format(#nbins,
                   lims,
                   lims,
                    "\""+xlab+"\"",
                    "\""+ylab+"\""))
        os.rename(r"rplot.png",name)
        
        if not save_file:
            del(name)
        os.remove("tmp.csv")    

        return Image(filename=name) 

    def compare_q_values(self, bm1,bm2, logx = False, logy = False, lims = "c(-1,10)", nbins = "30"):
        preds1 = bm1.get_q_val_table().reset_index()
        preds2 = bm2.get_q_val_table().reset_index()
        
        if len(preds1.columns)>2:
            preds1 = self.flatten_df(preds1)
        if len(preds2.columns)>2:
            preds2 = self.flatten_df(preds2)

        comparison_table = preds1.merge(preds2, on = "ProteinId")
        comparison_table.columns = ["id","x","y"]
        if logx:
            if comparison_table.x.max()<1:
                comparison_table.x = -1*np.log10(comparison_table.x)
        if logy:
            if comparison_table.y.max()<1:
                comparison_table.y = -1*np.log10(comparison_table.y)


        # if different signs
        if comparison_table.x.mean()*comparison_table.y.mean() < 0:
            comparison_table["x"] = -1*np.abs(comparison_table.x)
            comparison_table["y"] = -1*np.abs(comparison_table.y)

        return RPlotting().hexbin(comparison_table, "example_2.png","-log_10 adjusted P Value Perseus", "-log_10 adjusted P Value MD+ Discovery BYO", lims, nbins)
        
    def compare_log_FCs(self, bm1, bm2, nbins = 30):

        preds1 = bm1.get_estimated_logFCs().reset_index()
        preds2 = bm2.get_estimated_logFCs().reset_index()
        
        if len(preds1.columns)>2:
            preds1 = self.flatten_df(preds1)
        if len(preds2.columns)>2:
            preds2 = self.flatten_df(preds2)
     
        comparison_table = preds1.merge(preds2, on = "ProteinId")
        comparison_table.columns = ["id","x","y"]

        return RPlotting().hexbin(comparison_table, 
                                "example_1.png",
                                xlab = "log_2 fold change Perseus",
                                ylab ="log_2 fold change MD+ Discovery BYO",
                                lims = "c(-15,15)", 
                                nbins = nbins)
    
    def flatten_df(self,df):
        df.index = df.ProteinId
        df = df.drop("ProteinId", axis = 1)
        new_index = [str(i) + " " +str(c) for i in df.index for c in df.columns]
        values = df.values.flatten()
        new_df = pd.DataFrame(values,index = new_index)
        new_df["ProteinId"] = new_index
        new_df = new_df[["ProteinId",0]]
        return new_df