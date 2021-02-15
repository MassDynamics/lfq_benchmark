# Discovery Benchmark 

The discovery_benchmark codebase is designed to benchmark Label Free Quantification (LFQ) methods against benchmark datasets (iPRG2015) and the dynamic range benchmark dataset (referred to as UPS) dataset. 

For details about the methods enacted by this codebase please see: **our Discovery BYO Paper**

The codebase is made up of a few distinct components:

- benchmarkers: These classes organise information relating to LFQ output, theoretical and estimated predictions made about spiked proteins.
- ProteinTableLoader: This class is used to create a common protein table format amongst MD+ Discovery, MD+ Discovery BYO, Perseus and Proteome Discoverer Output. 
- Plotter: This class contains various plotting utilities. 
- RPlotting: This class uses Rpy2 module to run the hexbin package from R. 
- ConfusionMatrixCalculator: This class is used to generate confusion matrices using binary definitions of differential protein expression and detection. 
- compare_mqpars: This module enables users to quickly identify differences in maxquant parameter files. 

## Using this codebase:

Benchmarkers require the main protein quantification output to be contained inside a folder and the path to that folder to be specified when initializing the benchmarker. The software used to create this file must also be specified. 

For Example:
```
mq_home ="Path/to/Perseus/Output"
IPRG2015Benchmarker(mq_home, "Perseus").run()
```
## Example Output:

All of the figures contained in **BYO Paper** were produced using this package. Please see that document for examples. 

## Feedback:

This codebase was written mainly to benchark packages using numeric and visual representations of the correspondence between results and theoretical true values. 

As such it could be improved in many ways to make to easier to reuse/simpler/more elegant. 

If you have any ideas or suggestions, please share them and/or submit pull requests. 

## Versions and Requirements:

We used conda 4.9.2 and python 3.7.7 to manage packages required for this project. See the requirements.txt file for required packages (numpy, pandas, plotly, rpy2 and jupyter).

## Citations: 

Label Free Quantitation: Al Shweiki, M. R. et al. Assessment of Label-Free Quantification in Discovery Proteomics and Impact of Technological Factors and Natural Variability of Protein Abundance. J. Proteome Res. 16, 1410–1424 (2017).

iPRG2015: Choi, M. et al. ABRF Proteome Informatics Research Group (iPRG) 2015 Study: Detection of Differentially Abundant Proteins in Label-Free Quantitative LC-MS/MS Experiments. J. Proteome Res. 16, 945–957 (2017).

UPS: Cox, J. et al. Accurate proteome-wide label-free quantification by delayed normalization and maximal peptide ratio extraction, termed MaxLFQ. Mol. Cell. Proteomics 13, 2513–2526 (2014).

HER: Creedon, H. et al. Identification of novel pathways linking epithelial-to-mesenchymal transition with resistance to HER2-targeted therapy. Oncotarget 7, 11539–11552 (2016).

