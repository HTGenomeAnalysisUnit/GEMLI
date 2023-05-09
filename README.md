# GEMLI: Gene expression memory for lineage identification

GEMLI is an R package to predict cell lineages (cells with a common ancestor) from single cell RNA datasets (or any other type of single-cell gene expression data) and to call genes with a high gene expression memory on the predicted cell lineages. It is described in A.S. Eisele*, M. Tarbier*, A.A. Dormann, V. Pelechano, D.M. Suter | "Barcode-free prediction of cell lineages from scRNA-seq datasets" | 2022 bioRxiv: https://www.biorxiv.org/content/10.1101/2022.09.20.508646v1

The approach is based on findings of Phillips et al. 2019 (doi.org/10.1038/s41467-019-09189-8) where it was shown that some genes show varying gene expression across cell lineages that is stable over multiple cell generation.

## Installation
Simply run `library(devtools)` and then `install_github("UPSUTER/GEMLI", subdir="GEMLI_package_v0")`. If `devtools` is not installed yet, you can do so with `install.packages("devtools")`. GEMLI is now installed and can be used via `library(GEMLI)`. Additionally, GEMLI requires the packages `igraph` (visualization function) and `HiClimR` (fast correlation) loaded as `library(igraph)` and `library(HiClimR)`. To install these packages run `install.packages('igraph')` and `install.packages('HiClimR')` respectively.


## Development and feedback
We are still working to make GEMLI more intuitive, user-friendly, faster and versatile. Therefore existing functions will still be updated and new functionalities will be added. We'll publish a list of changes for each version for you to keep track. Your feedback is very welcome and will help us make GEMLI even better. What do you like about GEMLI? Something not working? What functions are you missing? Let us know! Contact: marcel.tarbier@scilifelab.se or almut.eisele@epfl.ch 

## Example 1: small lineages in mouse embryonic stem cells

For our first example we'll be looking at mouse embryonic stem cells that have been barcoded and cultured for 48h. We'll be working with a subset of this data for fast processing. In our subset we find 'family sizes' ranging from just two up to five related cells.

### Load example data
First we load the example data.

```
> load('GEMLI_example_data_matrix.RData')
> load('GEMLI_example_barcode_information.RData')
```

### Create a GEMLI items list
GEMLI's inputs and outputs are stored in a list of objects with predefined names. To run GEMLI you need at least a quality controlled and normalized gene expression matrix (rows = genes/features, colums = cells/samples). In this example we also provide a ground truth for lineages stemming from a barcoding experiment (values = barcode ID, names = cell IDs).

```
> GEMLI_items = list()
> GEMLI_items[['gene_expression']] = data_matrix
> GEMLI_items[['barcodes']] = lineage_dict_bc
>
> GEMLI_items[['gene_expression']][9:14,1:5]
                   AAACGAACAGGTGTGA-1 AAAGGTAGTTGCTTGA-1 AACCACAAGTTTGTCG-1 AAGCCATGTTCCACGG-1 AAGCGAGGTACGGCAA-1
ENSMUSG00000033845          14.761746         12.9570026          13.240645          8.8596794          12.791617
ENSMUSG00000025903           3.163231          1.2340002           4.878132          0.7383066           0.000000
ENSMUSG00000033813           6.326463          5.5530011           4.181256          8.8596794           5.482122
ENSMUSG00000002459           0.000000          0.0000000           0.000000          0.0000000           0.000000
ENSMUSG00000085623           0.000000          0.0000000           0.000000          0.0000000           0.000000
ENSMUSG00000033793           3.163231          0.6170001           1.393752          2.2149199           5.482122
>
> GEMLI_items[['barcodes']][1:5]
CACAGATAGTGATGGC-1 TATCTTGGTACGGGAT-1 AAACGAACAGGTGTGA-1 AGAGAATAGGTCATAA-1 GAGTGAGTCCAGTACA-1
                 2                  2                  2                  7                  7
```

### Perform lineage prediction
We can then identify cell lineages through repeated iterative clustering (this may take 2-3min). The `predict_lineages` function takes our GEMLI_items as input. It outputs a matrix of all cells against all cells with values corresponding to a confidence score that they are part of the same lineage. 

```
> GEMLI_items = predict_lineages(GEMLI_items)
>
> GEMLI_items[['prediction']][1:5,15:19]
                   AGAGAATAGGTCATAA-1 AGAGCAGCAAGTGATA-1 AGATGCTTCAAAGACA-1 AGGATCTGTATCGTTG-1 AGGGAGTAGACGATAT-1
AAACGAACAGGTGTGA-1                  0                  0                  0                  0                 14
AAAGGTAGTTGCTTGA-1                 87                  0                  0                  0                  0
AACCACAAGTTTGTCG-1                  0                  0                  0                  0                 39
AAGCCATGTTCCACGG-1                  0                  0                  0                  0                  0
AAGCGAGGTACGGCAA-1                  0                  0                  0                  0                  0
```

### Test lineage prediction
Since we have barcoding data for this dataset we can test the predicted lineages against our ground truth.
The `test_lineage_prediction` function again takes our `GEMLI_items` as input. It's important that a predcition has been run first with `predict_lineages`. It outputs the number of true positive predictions (TP), false positive predictions (FP), as well as precision and sensitivity for various confidence intervals. The output can be visualized by setting `plot_results` to `true`/`T`.

```
> GEMLI_items = test_lineages(GEMLI_items)
>
> GEMLI_items$testing_results
     TP    FP  precision sensitivity
0   274 24062 0.01125904   1.0000000
10  104   128 0.44827586   0.3795620
20   90    82 0.52325581   0.3284672
30   84    56 0.60000000   0.3065693
40   80    36 0.68965517   0.2919708
50   68    22 0.75555556   0.2481752
60   64    14 0.82051282   0.2335766
70   62    10 0.86111111   0.2262774
80   58     2 0.96666667   0.2116788
90   46     0 1.00000000   0.1678832
100  34     0 1.00000000   0.1240876
>
> test_lineages(GEMLI_items, plot_results=T)
```

<p align="center">
  <img width="500" height="500" src="https://github.com/UPSUTER/GEMLI/blob/main/Example/GEMIL_GitHub_testing.png">
</p>

### Visualize predictions as network
We can also investigate our predictions by visualizing them as a network with the `visualize_as_network` function. Here we need to set a `cutoff` that defines which predictions we want to consider. It represents a confidence score and high values yield fewer predictions with high precision while low values yield more predcitions with lower precision.

```
> visualize_as_network(GEMLI_items, cutoff=90) # left image
> visualize_as_network(GEMLI_items, cutoff=50) # right image
```
<p float="left">
  <img width="250" height="250" src="https://github.com/UPSUTER/GEMLI/blob/main/Example/GEMLI_GitHub_network_90.png">
  <img width="250" height="250" src="https://github.com/UPSUTER/GEMLI/blob/main/Example/GEMLI_GitHub_network_50.png">
</p>

If a ground truth e.g. from barcoding is avalable we can set `ground_truth` to `true`/`T` to highlight false predictions with red edges. Cells without barcode information will be displaye in white.

```
> visualize_as_network(GEMLI_items, cutoff=90, ground_truth=T) # left image
> visualize_as_network(GEMLI_items, cutoff=50, ground_truth=T) # right image
```
<p float="left">
  <img width="500" height="500" src="https://github.com/UPSUTER/GEMLI/blob/main/Example/GEMLI_GitHub_network_90_GT.png">
  <img width="500" height="500" src="https://github.com/UPSUTER/GEMLI/blob/main/Example/GEMLI_GitHub_network_50_GT.png">
</p>

### Extract lineage information
Now we can extract the lineage information with the `prediction_to_lineage_information` function. Again we need to set a `cutoff` that defines which predictions we want to consider. The function outputs both a lineage table and a 'dictionary', a vector that has the lineage number as values and the cell IDs as names.

```
> GEMLI_items = prediction_to_lineage_information(GEMLI_items, cutoff=50)
>
> GEMLI_items$predicted_lineage_table[1:5,]
     cell.ID              clone.ID
[1,] "AAACGAACAGGTGTGA-1" "1"
[2,] "AAAGGTAGTTGCTTGA-1" "2"
[3,] "AACCACAAGTTTGTCG-1" "3"
[4,] "AAGCCATGTTCCACGG-1" "4"
[5,] "AAGCGAGGTACGGCAA-1" "5"
>
> GEMLI_items$predicted_lineages[1:5,]
AAACGAACAGGTGTGA-1 AAAGGTAGTTGCTTGA-1 AACCACAAGTTTGTCG-1 AAGCCATGTTCCACGG-1 AAGCGAGGTACGGCAA-1
                 1                  2                  3                  4                  5
```

### Trim lineages that are too big

In some applications it may be useful to trim lineages that are too big. For instance if it is known that cells should have undergone only a certain number of divisions effectively limiting the lineage size or if you are only interested in sister cell pairs. Similarly, if you investigate large lineages you want to avoid lineages being merged due to few false predictions between otherwise well-interconnected lineages. The `suggest_network_trimming_to_size` function allows you to preview what a trimming to size would look like. It will again show the predicted lineages as networks but highlight all connections that would be trimmed given a certain size restriction (`max_size`). If you are happy with the suggested trimming you create new trimmed `GEMLI_items` list, in this example we called it `GEMLI_items_post_processed`. You can then again visualize the predictions to see how the changes have affected your predictions.
```
> suggest_network_trimming_to_size(GEMLI_items, max_size=2, cutoff=50) # left image
> GEMLI_items_post_processed = trim_network_to_size(GEMLI_items, max_size=2, cutoff=50)
> visualize_as_network(GEMLI_items_post_processed, cutoff=50) # right image
```
<p float="left">
  <img width="500" height="500" src="https://github.com/UPSUTER/GEMLI/blob/main/Example/GEMLI_GitHub_network_50_GT_ST.png">
  <img width="500" height="500" src="https://github.com/UPSUTER/GEMLI/blob/main/Example/GEMLI_GitHub_network_50_trim.png">
</p>


## Example 2: large lineages in intestinal crypts

In this example we'll be visualizing GEMLI results from a lineage-annotated dataset of murine intestinal crypts. Intestinal crypts originate from one or few intestinal stem cells, and are composed of a number of different cell types. In the data we find cells from individual crypts that range from just three up to forty-four related cells.

### Load example data.
First we load the example data. Here we already predicted the lineages using GEMLI and therefore do not include a count matrix, but rather start with the predictions right away.

```
> load('GEMLI_crypts_example_data_matrix.RData')
> load('GEMLI_crypts_example_barcode_information.RData')
```

### Create a GEMLI items list
We then create a GEMLI items list. This list is used to store the data, and create and store the outputs of GEMLI (for details check example one).

```
> GEMLI_items_crypts = list()
> GEMLI_items_crypts[['prediction']] = Crypts
> GEMLI_items_crypts[['barcodes']] = Crypts_bc_dict
```

### Visualize predictions as network
To visualize large lineages we'll use three different network layout algorithms: Fruchterman-Reingold, Kamada-Kawai, and grid. Each of them has advantages and disadvantages.

(1) Fruchterman-Reingold dispalys the cells of individual predicted organoids close together with ample space between them. This makes it hard to see connections within individual organoids but allows to get a good overview of individual structures.

(2) Kamada-Kawai spaces individual cells well, so we can see individual connections between them. It may, however, happen that two different predicted organoids are partially overlayed, as can be seen for dark red and bright green lineages on the right side of the plot.

(3) When the network is layed out as a grid, one gets generally a good overview of the predicted lineages and their connections, but it's hard to see which connections belongs to which cell in the same row.

```
> visualize_as_network(GEMLI_items_crypts, cutoff=70, display_orphan=F, max_edge_with=1, ground_truth=T, include_labels=F, layout_style="fr")
> visualize_as_network(GEMLI_items_crypts, cutoff=70, display_orphan=F, max_edge_with=1, ground_truth=T, include_labels=F, layout_style="kk")
> visualize_as_network(GEMLI_items_crypts, cutoff=70, display_orphan=F, max_edge_with=1, ground_truth=T, include_labels=F, layout_style="grid")
```

<p float="left">
  <img width="330" height="330" src="https://github.com/UPSUTER/GEMLI/blob/main/Example/GEMLI_GitHub_crypts_network_70_fr.png">
  <img width="330" height="330" src="https://github.com/UPSUTER/GEMLI/blob/main/Example/GEMLI_GitHub_crypts_network_70_kk.png">
  <img width="330" height="330" src="https://github.com/UPSUTER/GEMLI/blob/main/Example/GEMLI_GitHub_crypts_network_70_grid.png">
</p>

### Adding cell type information to the GEMLI items list
The cells of individual intestinal organoids can be assigned to different cell types. This information can be added to the GEMLI items list as 'cell_type' slot in the from of a dataframe with column 'cell.ID' and 'cell.type'.

```
> load('GEMLI_crypts_example_cell_type_annotation.RData')
> GEMLI_items[['cell_type']] = Crypts_annotation
```

### Color prediction network visualization by cell type
The visualization of the lineage predictions can now be colored by the cell type annotation. This allows to see the composition of individual intestinal organoids. 

```
> visualize_as_network(GEMLI_items, cutoff=70, max_edge_with=1, display_orphan=F, include_labels=F, ground_truth=T, highlight_FPs=T, layout_style="kk", cell_type_colors=T)

```
<p float="left">
  <img width="330" height="400" src="https://github.com/UPSUTER/GEMLI/blob/main/Example/GEMLI_GitHub_crypts_network_70_cell_type_colors.png">
</p>

Specific colors can be assigned to specific cell types by adding a dataframe with column 'cell.type' and 'color' in the GEMLI items list slot 'cell_type_color'.

```
> cell.type <- unique(GEMLI_items[['cell_type']]$cell.type)
> color <- c("#5386BD", "skyblue1", "darkgreen", "gold", "red", "darkred", "black")
> Cell_type_color <- data.frame(cell.type, color)
> GEMLI_items[['cell_type_color']] = Cell_type_color
>
> visualize_as_network(GEMLI_items, cutoff=70, max_edge_with=5, display_orphan=F, include_labels=F, ground_truth=T, highlight_FPs=T, layout_style="kk", cell_type_colors=T)

```

<p float="left">
  <img width="330" height="400" src="https://github.com/UPSUTER/GEMLI/blob/main/Example/GEMLI_GitHub_crypts_network_70_custom_cell_type_colors.png">
</p>




## Citation
If you use the package, please cite A.S. Eisele*, M. Tarbier*, A.A. Dormann, V. Pelechano, D.M. Suter | "Barcode-free prediction of cell lineages from scRNA-seq datasets" | 2022 bioRxiv: https://www.biorxiv.org/content/10.1101/2022.09.20.508646v1
