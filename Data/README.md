# Comparative benchmarking of single-cell clustering algorithms for transcriptomic and proteomic data

We employed 10 real datasets across 5 tissue types, encompassing over 50 cell types and more than 300,000 cells, each containing paired single-cell mRNA expression and surface protein expression data. These datasets were obtained using multi-omics technologies such as CITE-seq, ECCITE-seq, and Abseq. These paired multi-omics datasets, obtained by measuring gene or protein expression within the same set of cells, reflect identical biological conditions across the two omics. This consistency facilitates a comparable analysis of clustering algorithms, improving the reliability and comparability of the resulting analyses. The list of datasets used in this study is provided below.

## Datasets List
### 
|ID | Data                                                                 | Link                                                                 | Title                                                                                                               |PMID|
|:-------:|:-------:|:------------------------------------------------------------------------:|---------------------------------------------------------------------------------------------------------------------|:--:|
|1|Data1| [Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE143363) |Monocytic subclones confer resistance to venetoclax-based therapy in patients with acute myeloid leukemia.|[31974170](https://pubmed.ncbi.nlm.nih.gov/31974170/)|
|2|Data2| [Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148127) |Multi-modal Single-Cell Analysis Reveals Brain Immune Landscape Plasticity during Aging and Gut Microbiota Dysbiosis.  |[33264626](https://pubmed.ncbi.nlm.nih.gov/33264626/)|
|3|Data3| [Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163120) |Single-cell profiling of myeloid cells in glioblastoma across species and disease stage reveals macrophage competition and specialization.|[33782623](https://pubmed.ncbi.nlm.nih.gov/33782623/)|
|4|Data4| [Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163120) |Single-cell profiling of myeloid cells in glioblastoma across species and disease stage reveals macrophage competition and specialization.|[33782623](https://pubmed.ncbi.nlm.nih.gov/33782623/)|
|5|Data5| [Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128639) |Comprehensive integration of single-cell data. |[31178118](https://pubmed.ncbi.nlm.nih.gov/31178118/)|
|6|Data6| [Link](https://figshare.com/projects/Single-cell_proteo-genomic_reference_maps_of_the_human_hematopoietic_system/94469) |Single-cell proteo-genomic reference maps of the hematopoietic system enable the purification and massive profiling of precisely defined cell states.|[34811546](https://pubmed.ncbi.nlm.nih.gov/34811546/)|
|7|Data7| [Link](https://figshare.com/projects/Single-cell_proteo-genomic_reference_maps_of_the_human_hematopoietic_system/94469) |Single-cell proteo-genomic reference maps of the hematopoietic system enable the purification and massive profiling of precisely defined cell states.|[34811546](https://pubmed.ncbi.nlm.nih.gov/34811546/)|
|8|Data8| [Link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164378) |Integrated analysis of multimodal single-cell data. |[34062119](https://pubmed.ncbi.nlm.nih.gov/34062119/)|
|9|Data9| [Link1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE192742) [Link2](https://www.livercellatlas.org/index.php)   |Spatial proteogenomics reveals distinct and evolutionarily conserved hepatic macrophage niches.       |[35021063](https://pubmed.ncbi.nlm.nih.gov/35021063/)|
|10|Data10| [Link1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE192742) [Link2](https://www.livercellatlas.org/index.php) |Spatial proteogenomics reveals distinct and evolutionarily conserved hepatic macrophage niches.       |[35021063](https://pubmed.ncbi.nlm.nih.gov/35021063/)|

## DemoData
We provide a sample dataset for demonstration purposes, the `DemoData.rar` archive provides a sample dataset containing **1,000 cells**. 
After extracting the archive, the folder will contain the following files:

DemoData.rar  
├── ADT.csv  
├── ADTCount.csv  
├── celltypes.csv  
├── CITE.rds  
├── RNA.csv  
└── RNACount.csv  

To run the provided example script, place the extracted data into the designated dataset path as specified in the script. 

## Acknowledgments
We sincerely thank the researchers and institutions who provided the datasets used in this study. Their valuable contributions have made this work possible and continue to drive progress in the field.

