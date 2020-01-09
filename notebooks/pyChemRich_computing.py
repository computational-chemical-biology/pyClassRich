from pyChemRich import class_enrichment, impact_plot, ont_graph
import pandas as pd

chemrich = pd.read_csv('chemrich2classyfire/chemrich_2chemrich_input_file_with_clases.txt', sep='\t')

cls2rich = pd.read_csv('chemrich2classyfire/chemrich2classyfire.tsv', sep='\t')

chemrich_sup = chemrich.copy()
chemrich_sup['ChemRICHClusters'] = cls2rich['superclass_name']
chemrich_sub = chemrich.copy()
chemrich_sub['ChemRICHClusters'] = cls2rich['subclass_name']

clusterdf = class_enrichment(chemrich)
clusterdf_sub = class_enrichment(chemrich_sub)
clusterdf_sup = class_enrichment(chemrich_sup)
clusterdf_all = pd.concat([clusterdf_sub, clusterdf, clusterdf_sup]).reset_index(drop=True)

impact_plot(clusterdf, chemrich)
ont_graph(clusterdf_all)

gs_res = apply_gsea(taskid='3b1aa665-d10c-4712-b1d3-2b4d3d1659b5',
                    nap_taskid='5d9ddaf3be1b402bbf578e80348f0a2e',
                    mol_set='data/ChemOnt_2_1/dnp_sup_classes.json',
                    OUTDIR='data', method='ratio_of_classes')

