
import pandas as pd 
import numpy as np
import os 

import plotly.graph_objects as go
from plotly.graph_objects import Layout
import itertools
from itertools import combinations
from decimal import Decimal
import re

from sklearn.linear_model import LinearRegression
from scipy.stats import pearsonr

from datetime import datetime

class Plotter():

    def get_heatmap(self,df, fc = True):
        def df_to_plotly(df):
            return {'z': df.values.tolist(),
                    'x': df.columns.tolist(),
                    'y': df.index.tolist()}
        if fc:
            fig = go.Figure(data=go.Heatmap(df_to_plotly(df.sort_values(by = "ProteinId")), 
                            zmin = -4, zmax = 4, colorscale =[(0, "red"), (0.25, "yellow"), 
                                        (0.5, "green"), (0.75, "yellow"), (1, "red")],
                            ))
        else: 
            fig = go.Figure(data=go.Heatmap(df_to_plotly(df.sort_values(by = "ProteinId"))))
        print("True Positives:") 
        fig.show()
        self.export_fig(fig)
        return fig

    def scatter_log_fc_accuracy_by_protein(self, estLogFCs, trueLogFCs):
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=[-10, 10], y=[-10, 10], mode="lines",
                                    line=go.scatter.Line(
                                        color="gray", dash="dashdot"),
                                    showlegend=False))
        for i in estLogFCs.index:
            fig.add_trace(go.Scatter(x=list(trueLogFCs.loc[i]), y=list(estLogFCs.loc[i]), mode='markers',
                                    name=str(i), hovertext=str(i)))

        fig.update_layout(title= "LogFC Accuracy",
                        xaxis_title="True LogFC",
                        yaxis_title="Estimated LogFC")

        fig.show()
        self.export_fig(fig)
        return fig

    def scatter_log_fc_accuracy(self, estLogFCs, 
                                trueLogFCs, q_val_table):

        #convert to 1d and get labels
        ests = [list(estLogFCs.loc[i]) for i in estLogFCs.index]
        ests = list(itertools.chain.from_iterable(ests))
        trues = [list(trueLogFCs.loc[i]) for i in estLogFCs.index]
        trues = list(itertools.chain.from_iterable(trues))
        q_val = [list(q_val_table.loc[i]) for i in estLogFCs.index]
        q_val = list(itertools.chain.from_iterable(q_val))
        groups = [self.cf_group(q, est, pred) for (q, est, pred) in zip(q_val, ests, trues)]
        proteins = [list([i for x in estLogFCs.loc[i]]) for i in estLogFCs.index]
        proteins = list(itertools.chain.from_iterable(proteins))
        comparisons = [list(estLogFCs.loc[i].index) for i in estLogFCs.index]
        comparisons = list(itertools.chain.from_iterable(comparisons))
        labels = [p+" "+c+" "+"{:.2E}".format(Decimal(q)) for p,c,q in zip(proteins,comparisons,q_val)]

        data = pd.DataFrame({"ests":ests, "trues":trues, "qvals": q_val, "groups":groups, "proteins":proteins,
                            "comparisons":comparisons, "labels":labels})
        
        import plotly.graph_objects as go
        from plotly.graph_objects import Layout
        layout = Layout(plot_bgcolor='rgba(0,0,0,0)')
        fig = go.Figure(layout=layout)
        fig.add_trace(go.Scatter(x=[-15, 15], y=[-15, 15], mode="lines",
                                    line=go.scatter.Line(
                                        color="gray", dash="dashdot"),
                                    showlegend=False))
        
        colour_dict = {"TP":"Blue","FP":"Purple","TN":"Black", "FN":"Red"}
        for group in  ["TP", "FN", "TN", "FP"]:
            tmp = data[data.groups == group]
            fig.add_trace(go.Scatter(x=[None] + tmp.trues.to_list(), 
                                    y= [None] + tmp.ests.to_list(), 
                                    name=group,
                                    mode='markers',
                                    hovertext=labels,
                                    marker=dict(color=colour_dict[group]),
                                    #marker = dict(showscale = True,
                                    #colorbar=dict(title="p-value"), cmin = 0, cmax = 1),
                                    showlegend=True))

        reg = LinearRegression().fit(np.array(trues).reshape(-1, 1),np.array(ests))
        print("Intercept of a Linear Regression: ", reg.intercept_)
        print("Slope of a Linear Regression: ", reg.coef_)
        print("R2 of a Linear Regression: ", reg.score(np.array(trues).reshape(-1, 1),np.array(ests)))
        score = pearsonr(ests,trues)[0]
        title = "PearsonR = " + str(round(score,3))
        fig.update_layout(title= title,
                        xaxis_title="True Log2 Fold Change",
                        yaxis_title="Estimated Log2 Fold Change")

        fig.update_xaxes(showline=True, zeroline=True, 
                            linewidth=2, linecolor='black', 
                            gridcolor='Grey', mirror = True)
        fig.update_yaxes(showline=True, zeroline=True, 
                        linewidth=2, linecolor='black', 
                        gridcolor='Grey', mirror = True)

        fig.update_xaxes(zeroline=True, zerolinewidth=2, zerolinecolor='Grey')
        fig.update_yaxes(zeroline=True, zerolinewidth=2, zerolinecolor='Grey')

        fig.show()
        self.export_fig(fig)
        return fig

        import plotly.graph_objects as go

    def volcano_plot(self, protein_table, comparison = None, logy = False):
        
        if comparison:
            protein_table = protein_table.filter(regex = comparison)
            protein_table = protein_table.rename(columns=lambda x: re.sub(comparison,'',x))
        
        #create figure
        fig = go.Figure()
        x = protein_table["logFC"]
        y = -1*np.log10(protein_table["adj.P.Val"])

        width = max(np.abs(x))
        fig.add_trace(go.Scatter(x=x, y=y,marker=dict(size=5,
                                    line=dict(width=1,
                                                color='DarkSlateGrey')),
                                mode='markers',
                                hovertext=protein_table.index, showlegend=False))

        fig = self.format_volcano_plot(fig, width)
        fig.show()
        self.export_fig(fig)
        return fig

    def format_volcano_plot(self, fig, width):
        fig.add_trace(go.Scatter(x=[-1*width, width], y=[-1*np.log10(0.05), -1*np.log10(0.05)], mode="lines",
                                    line=go.scatter.Line(
                                        color="gray", width = 1),
                                    showlegend=False))
        fig.add_trace(go.Scatter(x=[-1, -1], y=[0,-1*np.log10(10**(-13))], mode="lines",
                                    line=go.scatter.Line(
                                        color="gray", dash="dash"),
                                    showlegend=False))
        fig.add_trace(go.Scatter(x=[1, 1], y=[0,-1*np.log10(10**(-13))], mode="lines",
                                    line=go.scatter.Line(
                                        color="gray", dash="dash"),
                                    showlegend=False))

        fig.add_trace(go.Scatter(x=[width-3],
                                y=[-1*np.log10(0.05)+0.5],
                                mode="text",
                                name="Markers and Text",
                                text=["5% FDR Line"],
                                textposition="bottom center", 
                                showlegend=False))

        fig.update_xaxes(showgrid=False)
        fig.update_yaxes(showgrid=False)
        #fig.update_layout(plot_bgcolor='rgba(0,0,0,0)')

        fig.update_layout(
            xaxis_title = "Log2 Fold Changes",
            yaxis_title = "-log10 Adjusted P-value",
            xaxis = dict(
                tickmode = 'array',
                tickvals = [-10, -8, -6, -4, -2, 0, 2,4,6,8,10],
                ticktext =  [-10, -8, -6, -4, -2, 0, 2,4,6,8,10]
            )
        )


        return fig

    def compare_log_FCs(self, bm1, bm2):

        preds1 = bm1.get_estimated_logFCs().reset_index()
        preds2 = bm2.get_estimated_logFCs().reset_index()

        if len(preds1.columns)>2:
            preds1 =  Plotter().flatten_df(preds1)
        if len(preds2.columns)>2:
            preds2 =  Plotter().flatten_df(preds2)

        comparison_table = preds1.merge(preds2, on = "ProteinId")

        layout = Layout(plot_bgcolor='rgba(0,0,0,0)')
        fig = go.Figure(layout = layout)
        fig.add_trace(go.Scatter(x=[-15, 15], y=[-15, 15], mode="lines",
                                    line=go.scatter.Line(
                                        color="gray", dash="dashdot"),
                                    showlegend=False))

        fig.add_trace(go.Scatter(x=comparison_table.iloc[:,1], 
                                y=comparison_table.iloc[:,2], 
                                marker=dict(color="black"),
                                mode='markers',
                                hovertext=comparison_table.iloc[:,0],
                                showlegend=False))
                
        score = pearsonr(comparison_table.iloc[:,1],comparison_table.iloc[:,2])[0]
        title = "PearsonR = " + str(round(score,3))
        fig.update_layout(title= title,
                        xaxis_title="Perseus: Estimated Fold Change",
                        yaxis_title="MD 1.0: Estimated Fold Change")

        fig.update_xaxes(showline=True, zeroline=True, 
                            linewidth=2, linecolor='black', 
                            gridcolor='Grey', mirror = True)
        fig.update_yaxes(showline=True, zeroline=True, 
                        linewidth=2, linecolor='black', 
                        gridcolor='Grey', mirror = True)

        fig.update_xaxes(zeroline=True, zerolinewidth=2, zerolinecolor='Grey')
        fig.update_yaxes(zeroline=True, zerolinewidth=2, zerolinecolor='Grey')

        fig.write_html("logFCComparison.html")
        self.export_fig(fig)
        return fig

    def compare_q_values(self, bm1, bm2):

        preds1 = bm1.get_q_val_table().reset_index()
        preds2 = bm2.get_q_val_table().reset_index()
        
        if len(preds1.columns)>2:
            preds1 = self.flatten_df(preds1)
        if len(preds2.columns)>2:
            preds2 = self.flatten_df(preds2)
        
        comparison_table = preds1.merge(preds2, on = "ProteinId")
        comparison_table["Proteins"] = comparison_table.ProteinId.apply(lambda x: x.split(" ")[0])



        layout = Layout(plot_bgcolor='rgba(0,0,0,0)')
        fig = go.Figure(layout = layout)
        fig.add_trace(go.Scatter(x=[0, 10], y=[0, 10], mode="lines",
                                    line=go.scatter.Line(
                                        color="gray", dash="dashdot"),
                                    showlegend=False))

        fig.add_trace(go.Scatter(x=-1*np.log10(comparison_table.iloc[:,1]), 
                                y=-1*np.log10(comparison_table.iloc[:,2]), 
                                marker=dict(color="black"),
                                mode='markers',
                                hovertext=comparison_table.iloc[:,0],
                                showlegend=False))


        score = pearsonr(comparison_table.iloc[:,1],comparison_table.iloc[:,2])[0]
        title = "PearsonR = " + str(round(score,3))
        fig.update_layout(title= title,
                        xaxis_title="Perseus: -log10(adjusted P Value)",
                        yaxis_title="MD 1.0: -log10(adjusted P Value)")

        fig.update_xaxes(showline=True, zeroline=True, 
                            linewidth=2, linecolor='black', 
                            gridcolor='Grey', mirror = True)
        fig.update_yaxes(showline=True, zeroline=True, 
                        linewidth=2, linecolor='black', 
                        gridcolor='Grey', mirror = True)

        fig.update_xaxes(zeroline=True, zerolinewidth=2, zerolinecolor='Grey')
        fig.update_yaxes(zeroline=True, zerolinewidth=2, zerolinecolor='Grey')

        fig.write_html("QValueComparison.html")
        self.export_fig(fig)
        return fig

    def flatten_df(self,df):
        df.index = df.ProteinId
        df = df.drop("ProteinId", axis = 1)
        new_index = [str(i) + " " +str(c) for i in df.index for c in df.columns]
        values = df.values.flatten()
        new_df = pd.DataFrame(values,index = new_index)
        new_df["ProteinId"] = new_index
        new_df = new_df[["ProteinId",0]]
        return new_df

    def cf_group(self, q_val, estimate, true_value):

        difference =  np.abs(true_value) >= 1
        same_sign_detection = (estimate>0)==(true_value>0) and (q_val <0.05)

        if difference and same_sign_detection:
            return "TP"
        elif ~difference and same_sign_detection:
            return "FP"
        elif difference and ~same_sign_detection:
            return "FN"
        elif ~difference and ~same_sign_detection:
            return "TN"

    def export_fig(self, fig):
        if not os.path.exists("plots"):
            os.mkdir("plots")

        now = datetime.now()
        current_time = now.strftime("%H_%M_%S")
        fig.write_image("plots/fig1"+current_time+".png", engine="kaleido", scale = 2)
        return 