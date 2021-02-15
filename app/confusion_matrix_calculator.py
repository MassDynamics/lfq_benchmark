from sklearn.metrics import confusion_matrix
import pandas as pd

class ConfusionMatrixCalculator():

    def run(self, bm):

        comparison = self.get_comparison_table(bm)

        return confusion_matrix(comparison.Real,comparison.Positive)

    def get_differentially_expressed_proteins(self, theoretical_log_FCs):

        cols = [col for col in theoretical_log_FCs.columns if "logFC" in col]
        dfs = []
        for col in cols:

            comparison = col[5:]
            #copy table
            df = theoretical_log_FCs.copy()

            cols_regex = col
            df = df.filter(regex = cols_regex , axis = 1)

            #rename index
            df.index = df.index +" " + comparison
            
            #rename columns
            df.columns = ["logFC"]
            
            dfs.append(df)

        reals = pd.concat(dfs)
        reals["Real"] = ((reals.logFC > 1) | (reals.logFC < -1))
        return reals

    def get_detections(self,protein_table):

        #for key, val in spike_mapping.items():
        #    protein_table.index = protein_table.index.str.replace(val,key)
        cols = [col for col in protein_table.columns if "logFC" in col]
        dfs = []
        for col in cols:

            comparison = col[5:]
            #copy table
            df = protein_table.copy()

            cols_regex = comparison
            df = df.filter(regex = cols_regex , axis = 1)

            #rename index
            df.index = df.index +" " + comparison
            
            #rename columns
            df = df.reindex(sorted(df.columns), axis=1)
            df.columns = ["adj.P.Val", "logFC"]
            
            dfs.append(df)

        detections = pd.concat(dfs)

        #print("Max p value: ", detections["adj.P.Val"].max())
        if detections["adj.P.Val"].max() > 1:
            detections["adj.P.Val"] = 10**(-1*detections["adj.P.Val"])

        detections["Positive"] = (detections["adj.P.Val"]<0.05) & ((detections.logFC > 1) | (detections.logFC < -1))

        return detections
    
    def get_comparison_table(self, bm):

        protein_table = bm.protein_table
        theoretical_log_FCs = bm.get_theoretical_logFCs()

        det_expressed_proteins = self.get_detections(protein_table)
        dif_expressed_proteins = self.get_differentially_expressed_proteins(theoretical_log_FCs)
        dif_expressed_proteins.index = dif_expressed_proteins.index.str.replace("UPS", "")
        det_expressed_proteins.index = det_expressed_proteins.index.str.replace("UPS", "")
        det_expressed_proteins.index = det_expressed_proteins.index.str.replace("sampleM","1")
        det_expressed_proteins.index = det_expressed_proteins.index.str.replace("sampleL","2")
        det_expressed_proteins.index = det_expressed_proteins.index.str.replace("sample","")
        dif_expressed_proteins.index = dif_expressed_proteins.index.str.replace("sample","")
     
        comparison = pd.merge(det_expressed_proteins,dif_expressed_proteins, how = "left", on = "ProteinId")
        comparison.index = det_expressed_proteins.index
        comparison.Real = comparison.Real.fillna(False)

        for i in dif_expressed_proteins.index:
            if i not in comparison.index:
                print("comparison {} wasn't matched".format(i))

        return comparison