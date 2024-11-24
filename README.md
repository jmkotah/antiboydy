# antiboydy
R package of functions used during my PhD

Last updated: 2024-11-24

Required packages: dplyr, ggplot2, ggrepel

To install:
1) Make sure you have the devtools package installed in R: <i>install.packages("devtools")</i>
2) Install this package with the following command: <i>devtools::install_github("https://github.com/jmkotah/antiboydy/")</i>
3) Load package <i>library(antiboydy)</i>

For suggestions/comments, contact me at j.m.kotah@umcg.nl or open an issue.

Currently available functions:
1) JK_outlier: detects outliers using the IQR method
2) ES_bw_gain: used to calculate P9-P2 body weight in early-life stress experiments, using average P2 weights per nest per sex in place of individual P2 weights
3) jk_volcano: create volcano plots from any dataset with foldchange and p-val/q-val/FDR information
4) jk_sharedDE and jk_4wayplot: input two differential expression dataframes and make four way plots