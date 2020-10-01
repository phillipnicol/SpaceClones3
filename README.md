# SpaceClones3

This is an R package which implements a statistical procedure for finding tumor subclones from spatially resolved RNA-seq data (i.e., spatial transcriptomics, Slide-seq, FISH).
This method is described in a conference paper (which will be available soon). This is currently a very minimal 
R package containing a minimal implementation of the algorithm. If there is interest, I would be happy to improve the code and documentation. 

## Installation 
Open R console and run `devtools::install_github("phillipnicol/SpaceClones3")`. 

## Input 

All that is needed a gene expression data frame `X` and a spatial weights matrix `D`. The rows of `X` corresponds to sites and the columns are genes. Entries in `X` should be the log2 fold change
of that gene compared to the mean expression of the control sample. `X` should only contain cells/sites which are believed to be cancer. The square matrix `D` gives the spatial weights between cells. In the paper I suggest using the inverse of the
squared Euclidean distance between the two cells, but symmetric measure can be used here. The remaining optional arguments pertain to cutoffs described in the paper (none of these will make sense if you haven't read the paper, just use the default values). 

## Output 

The output is an object of class [Mercator](https://CRAN.R-project.org/package=Mercator). Of particular interest will be `out$clusters` which gives an assignment of genes to subclones. 
