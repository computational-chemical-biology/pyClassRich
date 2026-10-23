# pyClassRich 

<p align="center">
  <img src="https://github.com/computational-chemical-biology/pyClassRich/blob/master/img/pyclassrich.jpeg" alt="logo"/>
</p>

pyClassRich uses class annotation on LC-MS experiments to perform enrichment analysis, generating ontology weighted graphs and impact plots.

### Install conda

```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```
   
### Create a dedicated conda environment and activate

```
conda env create -f environment.yml
conda activate pyclassrich
pip install git+https://github.com/computational-chemical-biology/pyClassRich.git
```

### Alternative installtion 

```
conda env create -f environment_alt.yml
conda activate pyclassrich
```

Download the package and inside the package directory do

```
python setup.py install

```


### License

This project is licensed under the BSD 3-Clause License - see the [LICENSE](LICENSE) file for details

