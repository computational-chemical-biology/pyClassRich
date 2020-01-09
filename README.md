# lprojection 
Projection layout for graphs.

# Installation
```
conda create -n lproj python=3.5   
conda install -n lproj rdkit -c rdkit 
source activate lproj
pip install git+https://gitlab.com/rsilvabioinfo/lprojection.git
```
# Getting started

Run the following to access the help function

```
layout_script layout --help
```

# Creating class colored tSNE layout 

Run the following command to create a projection based layout gml file, colored by compound class. 

```
layout_script layout \
    --taskid '4fd90fa52ade4bdcb604945952964fa0' \
    --projection TSNE \
    --workflow V2 \
    --scaling 30 \
    --linput Cosine \
    --cthr 0.7 \
    --meta data/associated_meta_inchi_match_full.tsv \
    --metac superclass_name

```
# Parameter optimization 

Run the following command to optimize the tSNE parameters. 

```
layout_script tsne-optim \
    --taskid '4fd90fa52ade4bdcb604945952964fa0,608e5eb3402f4dd599fcc88b1c8a40e9' \
    --projection TSNE \
    --workflow V2 \
    --meta data/associated_meta_inchi_match_full.tsv 
```

