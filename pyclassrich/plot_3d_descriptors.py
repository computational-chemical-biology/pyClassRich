import pandas as pd
import numpy as np
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors

import skbio
from qiime2.plugins import metadata, feature_table, diversity, emperor
from q2_emperor import plot, procrustes_plot, biplot, generic_plot
import qiime2
from qiime2 import Artifact

def smiles2inchi(smiles):
    return Chem.MolToInchi(Chem.MolFromSmiles(smiles))

def filter_extremes(data, lower_quant,
                    higher_quant):
    pc1 = (data['PC1'] > data['PC1'].quantile([lower_quant]).values[0]) & (data['PC1'] < data['PC1'].quantile([higher_quant]).values[0])
    pc2 = (data['PC2'] > data['PC2'].quantile([lower_quant]).values[0]) & (data['PC2'] < data['PC2'].quantile([higher_quant]).values[0])
    pc3 = (data['PC3'] > data['PC3'].quantile([lower_quant]).values[0]) & (data['PC3'] < data['PC3'].quantile([higher_quant]).values[0])
    return data[pc1 & pc2 & pc3]

def get_descr(data, inchi_col, smiles_col=None,
              nms = ['MolLogP', 'ExactMolWt', 'MolMR']):
    if smiles_col!=None:
        inchi = []
        for i in data.index:
            if data.loc[i, smiles_col]!='' and type(data.loc[i, inchi_col])==float:
                try:
                    inchi.append(smiles2inchi(data.loc[i, smiles_col]))
                except:
                    inchi.append('')

        data.loc[(data.loc[i, smiles_col]!='') & (data[inchi_col].isna()), inchi_col] = inchi

    mols = [Chem.MolFromInchi(x) for x in data[inchi_col]]
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(nms)
    descrs = [calc.CalcDescriptors(x) for x in mols]
    descrs = pd.DataFrame(descrs, columns=nms)
    return descrs


# See Descriptors._descList
def create_emperor_plot(data, meta, data_index, meta_index,
                        out_dir, scale=False):

    if data_index in data.columns:
        data.set_index(data_index, inplace=True)

    if meta_index in meta.columns:
        meta.set_index(meta_index, inplace=True)
    meta.index.name = 'id'
    meta = meta[~meta.index.isnull()]

    data = data[data.index.isin(meta.index)]
    if scale:
        data = data.apply(lambda a: (a-a.mean())/a.std())

    eigvals = pd.Series(np.array([0.50, 0.25, 0.25]),
                        index=['PC1', 'PC2', 'PC3'])
    proportion_explained = pd.Series([0.34, 0.33, 0.33], index=['PC1', 'PC2', 'PC3'])
    data.columns = ['PC1', 'PC2', 'PC3']

    #data['PC1'].quantile([0.05, 0.1, 0.25, 0.5, 0.75])

    pcoa = skbio.OrdinationResults('PCoA',
           'Principal Coordinate Analysis',
            eigvals,
            data,
            proportion_explained=proportion_explained)

    metadata = qiime2.Metadata(meta)
    pcoa = Artifact.import_data('PCoAResults', pcoa)
    emperor_plot = emperor.visualizers.plot(pcoa,
                                            metadata)
    if '.qzv' in out_dir:
        emperor_plot.visualization.save(out_dir)
    else:
        emperor_plot.visualization.export_data(out_dir)

