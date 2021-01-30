#https://github.com/qiime2/q2-emperor/blob/master/q2_emperor/tests/test_plot.py
import qiime2
from qiime2 import Artifact
import pandas as pd
from qiime2.plugins import metadata, feature_table, diversity, emperor
from scipy.spatial.distance import squareform, pdist
import skbio
from q2_emperor import plot, procrustes_plot, biplot, generic_plot

dm1 = pd.read_csv('iqt2_mldist.tsv', sep='\t')
dm2 = pd.read_csv('euclidean_scaled_metabolomis.tsv', sep='\t')
meta = pd.read_csv('Tabela_IdentificacaoCodigos_36Pseudonocardia.csv', sep='\t')

dm1 = skbio.DistanceMatrix(dm1.to_numpy(), ids=map(str, dm1.index.tolist()))
dm1 = Artifact.import_data("DistanceMatrix", dm1)
pcoa1 = diversity.methods.pcoa(dm1)

dm2 = skbio.DistanceMatrix(dm2.to_numpy(), ids=map(str, dm2.index.tolist()))
dm2 = Artifact.import_data("DistanceMatrix", dm2)
pcoa2 = diversity.methods.pcoa(dm2)

prcst = diversity.methods.procrustes_analysis(pcoa1.pcoa, pcoa2.pcoa)

meta.drop(meta.columns[0], axis=1, inplace=True)
meta.rename(columns={meta.columns[0]:'ID'}, inplace=True)
meta.rename(columns={meta.columns[1]:'Name'}, inplace=True)
meta.drop([36], inplace=True)
meta.drop([6], inplace=True)
meta.reset_index(drop=True, inplace=True)
meta['Name'] = meta['Name'].str.replace(' ICBG', '')
meta.set_index('ID', inplace=True)
metadata = qiime2.Metadata(meta)

emperor_plot = emperor.visualizers.plot(pcoa1.pcoa,
                                        metadata)
emperor_plot.visualization.export_data('phylo_pcoa')
emperor_plot.visualization.save('phylo_pcoa.qzv')

emperor_plot = emperor.visualizers.plot(pcoa2.pcoa,
                                        metadata)
emperor_plot.visualization.export_data('metabo_pcoa')
emperor_plot.visualization.save('metabo_pcoa.qzv')

procrustes_plot = emperor.visualizers.procrustes_plot(pcoa1.pcoa,
                                                     other_pcoa=pcoa2.pcoa,
                                                     metadata=metadata)
procrustes_plot.visualization.export_data('procrustes')
procrustes_plot.visualization.save('procrustes.qzv')
