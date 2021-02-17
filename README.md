## Installation
<p align="center">
  <img src="https://gitlab.com/rsilvabioinfo/chemwalker/-/blob/master/img/walker.gif" alt="logo"/>
</p>

pyClassRich uses class annotation on LC-MS experiments to perform enrichment analysis, generating ontology weighted graphs and [GSEA](https://github.com/ostrokach/gseapy) results.

### Install conda

```
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh

```
   
### Create a dedicated conda environment and activate

```
conda env create -f environment.yml
pip install git+https://gitlab.com/rsilvabioinfo/pyclassrich.git
```
