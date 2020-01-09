
import pandas as pd
import numpy as np
from collections import Counter
from scipy import stats
from statsmodels.stats.multitest import multipletests
import gseapy as gp

import requests
import io
import yaml
import os
import json
import xmltodict
from rdkit import Chem

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import matplotlib.colors as mcolors
import matplotlib as mpl
plt.style.use('ggplot')
#mpl.rcParams['lines.markersize'] = 200

import networkx as nx
import obonet

def apply_gsea(taskid, nap_taskid, mol_set, OUTDIR, method='ratio_of_classes'):
    """ Applies GSEA to contrast defined groups on FBMN feature tables
    Parameters
    ----------
    taskid : str
       Carcara task id
    nap_taskid : str
       NAP task id
    mol_set : str or dict
       Molecular class dictionary
    method : str
      GSEA methods: 'signal_to_noise', 't_test', 'ratio_of_classes',
      'diff_of_classes', 'log2_ratio_of_classes'
    Returns
    -------
    gseapy : object
       gseapy object or dict for multiple comparisons
    """
    if type(mol_set)==str:
        with open(mol_set, 'r') as f:
            mol_set = json.load(f)

    dout = os.path.join(OUTDIR, taskid)
    with open(os.path.join(dout, 'params.yml'), 'r') as stream:
        params = yaml.safe_load(stream)

    url_to_param = 'http://proteomics2.ucsd.edu/ProteoSAFe/DownloadResultFile?task=%s&block=main&file=params/' % nap_taskid
    r = requests.get(url_to_param)
    plist = xmltodict.parse(r.text)['parameters']['parameter']
    gtaskid = [x['#text'] for x in plist if x['@name']=='JOBID'][0]

    url_to_metadata = "http://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task=%s&block=main&file=metadata_table/" % gtaskid
    url_to_quantification = "http://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task=%s&block=main&file=quantification_table/" % gtaskid

    annotated = pd.read_csv(os.path.join(dout,'annotated_gnps_nap.tsv'),
                            sep='\t')

    meta = pd.read_csv(io.StringIO(requests.get(url_to_metadata).text), sep='\t')
    feat = pd.read_csv(io.StringIO(requests.get(url_to_quantification).text))

    feat = pd.merge(annotated[['cluster.index', 'InchiKey']], feat, left_on='cluster.index', right_on='row ID')
    feat.drop(['cluster.index', 'row m/z', 'row retention time'], axis=1, inplace=True)
    # check if all inchikeys are calculated
    feat.rename(columns={'InchiKey':'NAME', 'row ID':'DESCRIPTION'}, inplace=True)
    feat.columns = feat.columns.str.replace('(.+)\\.mzX?ML .+', '\\1')
    meta['filename'] = meta['filename'].str.replace("\\.mzXML|\\.mzML", "")

    if 'filter' in params.keys():
        if params['filter']['perform']:
            fvalues = params['filter']['value'].split(',')
            meta = meta[~meta[params['filter']['field']].isin(fvalues)]
    if 'keep' in params.keys():
        if params['keep']['perform']:
            fvalues = params['keep']['value'].split(',')
            meta = meta[meta[params['keep']['field']].isin(fvalues)]

    cn = params['comparison']['field']
    if params['comparison']['classes']=='all':
        cs = meta[cn].value_counts()
        cs = cs[cs>2].index.tolist()
    else:
        cs = params['comparison']['classes'].split(',')

    meta = meta[meta[cn].isin(cs)]
    cls_vec = meta[cn].tolist()
    feat = feat[['NAME', 'DESCRIPTION']+meta['filename'].tolist()]
    feat = feat[~((feat['NAME']=='') | (feat['NAME'].isnull()))]
    feat['DESCRIPTION'] = 0

    ex = feat[feat.columns[2:]].apply(lambda a: sum(a)==0, axis=1)
    if ex.sum():
        feat = feat[~ex]

    if params['pcoa_norm']:
        feat[feat.columns[2:]] = feat[feat.columns[2:]].apply(lambda a: a/sum(a))

    #if params['pcoa_scale']:
    #    df2 = (df2-df2.mean())/df2.std()
    # Average duplicates, how can we reduce duplicates?
    feat = feat.groupby('NAME', as_index=False).mean()
    feat['DESCRIPTION'] = 'na'

    if len(cs)>2:
        gsea_dict = {}
        for c in cs:
            cls_vec_tmp = np.array(cls_vec)
            cls_vec_tmp[cls_vec_tmp!=c] = 'others'
            gs_res = gp.gsea(data=feat, # or data='./P53_resampling_data.txt'
                             gene_sets=mol_set, # enrichr library names
                             cls=cls_vec_tmp.tolist(), # './data/P53.cls'
                             # set permutation_type to phenotype if samples >=15
                             permutation_type='phenotype',
                             permutation_num=1000, # reduce number to speed up test
                             outdir=None,  # do not write output to disk
                             no_plot=True, # Skip plotting
                             method=method,
                             processes=4,
                             format='png')
            gsea_dict[c] = gs_res
        return gsea_dict
    else:
        gs_res = gp.gsea(data=feat, # or data='./P53_resampling_data.txt'
                         gene_sets=mol_set, # enrichr library names
                         cls=cls_vec, # './data/P53.cls'
                         # set permutation_type to phenotype if samples >=15
                         permutation_type='phenotype',
                         permutation_num=1000, # reduce number to speed up test
                         outdir=None,  # do not write output to disk
                         no_plot=True, # Skip plotting
                         method=method,
                         processes=4,
                         format='png')
        return gs_res

def kstest(df, clusterid, cfield='ChemRICHClusters', pfield='pvalue',
          alpha=0.05):
    clmember = np.where(df[cfield]==clusterid)[0]
    pvalmembers = df.loc[clmember, pfield]
    if sum(pvalmembers<alpha) > 1:
        pres = stats.kstest(pvalmembers, 'uniform', alternative = 'greater')
        return pres[1]
    else:
        return 1


#clusterids = chemrich['ChemRICHClusters'].unique().tolist()

def class_enrichment(chemrich, cfield='ChemRICHClusters', nfield='CompoundName', pfield='pvalue'):
    chemdict = Counter(chemrich[cfield])
    clusterids = np.array(list(chemdict.keys()))
    clusterids = clusterids[np.array(list(chemdict.values()))>2]
    clusterids = clusterids[(clusterids!='nan') & (clusterids!='')]

    clspvals = [kstest(chemrich, x, cfield=cfield, pfield=pfield) for x in clusterids]
    clspvals = np.array(clspvals)
    # May be not needed
    clspvals[clspvals==0] = 2**(-20)
    clusterdf = pd.DataFrame([clusterids, clspvals], index=['name', 'pvalues']).T

    keycpdname = []
    for x in clusterdf['name']:
        keycpdname.append(chemrich.loc[chemrich.loc[chemrich[cfield]==x, 'pvalue'].idxmin(), nfield])

    clusterdf['keycpdname'] = keycpdname

    altrat = []
    for x in clusterdf['name']:
        alt = sum((chemrich[cfield]==x) & (chemrich['pvalue'] < 0.05))/sum(chemrich[cfield]==x)
        altrat.append(alt)

    uprat = []
    for x in clusterdf['name']:
        up = sum((chemrich[cfield]==x) & (chemrich['pvalue'] < 0.05) & (chemrich['foldchange'] > 1))
        upsum = sum((chemrich[cfield]==x) & (chemrich['pvalue'] < 0.05))
        if upsum:
            up = up/upsum
        else:
            up = np.nan
        uprat.append(up)

    clustsvec = [sum(chemrich[cfield]==x) for x in clusterdf['name']]

    alteredMetabolites = []
    for x in clusterdf['name']:
        alt = sum((chemrich[cfield]==x) & (chemrich['pvalue'] < 0.05))
        alteredMetabolites.append(alt)

    clusterdf['alteredMetabolites'] = alteredMetabolites

    upcount = []
    for x in clusterdf['name']:
        up = sum((chemrich[cfield]==x) & (chemrich['pvalue'] < 0.05) & (chemrich['foldchange'] > 1))
        upcount.append(up)

    clusterdf['upcount'] = upcount

    downcount = []
    for x in clusterdf['name']:
        down = sum((chemrich[cfield]==x) & (chemrich['pvalue'] < 0.05) & (chemrich['foldchange'] < 1))
        downcount.append(down)

    clusterdf['downcount'] = downcount

    clusterdf['upratio'] = uprat
    clusterdf['altratio'] = altrat
    clusterdf['csize'] = clustsvec

    clusterdf = clusterdf[clusterdf['csize'] > 2]
    clusterdf['adjustedpvalue'] = multipletests(clusterdf['pvalues'], method = "fdr_bh")[1]

    xlogp = []
    for x in clusterdf['name']:
        logp = np.median(chemrich.loc[chemrich[cfield]==x, 'xlogp'])
        xlogp.append(logp)

    clusterdf['xlogp'] = xlogp
    return clusterdf


def label_point(x, y, val, ax):
    a = pd.concat({'x': x, 'y': y, 'val': val}, axis=1)
    for i, point in a.iterrows():
        ax.text(point['x']+.02, point['y'], str(point['val']), fontsize=9)


def impact_plot(clusterdf, chemrich):
    Compounds = []
    for x in clusterdf['name']:
        cpds = chemrich.loc[chemrich['ChemRICHClusters']==x, 'CompoundName']
        Compounds.append('<br>'.join(cpds.tolist()))

    clusterdf['Compounds'] = Compounds

    clustdf = clusterdf[clusterdf['pvalues'] != 1 ]
    clustdf_alt_impact = clustdf[(clustdf['pvalues'] < 0.05) &
                                  (clustdf['csize'] > 1) &
                                  (clustdf['alteredMetabolites'] > 1)]

    clustdf_alt_impact.sort_values(['xlogp'], inplace=True)
    clustdf_alt_impact['order'] = range(clustdf_alt_impact.shape[0])
    clustdf_alt_impact['logPval'] = [-np.log(x) for x in clustdf_alt_impact['pvalues']]

    cmap = {obsv_id:clustdf_alt_impact.loc[obsv_id, "upratio"] for obsv_id in clustdf_alt_impact.index}

    clist = list(sns.color_palette("RdBu_r", 7).as_hex())
    #clist.reverse()
    sm = ScalarMappable(norm=Normalize(vmin=min(list(cmap.values())),
                                       vmax=max(list(cmap.values()))),
                        cmap=ListedColormap(clist))


    s = [20*n for n in clustdf_alt_impact["logPval"]]

    fig, ax = mpl.pyplot.subplots()
    xlogp = clustdf_alt_impact["xlogp"]
    logPval = clustdf_alt_impact["logPval"]
    ax.scatter(x="xlogp", y="logPval",
               s=s,
               color=[sm.to_rgba(cmap[obsv_id]) for obsv_id in
                      clustdf_alt_impact.index],
               data=clustdf_alt_impact)

    label_point(clustdf_alt_impact["xlogp"], clustdf_alt_impact["logPval"],
                clustdf_alt_impact["name"],
                plt.gca())

    xmin, xmax, ymin, ymax = ax.axis()
    ax.axis([min(xlogp.min()-0.3*xlogp.min(), xmin), xlogp.max()+0.3*xlogp.max(),
              min([logPval.min()-0.3*logPval.min(), ymin]), logPval.max()+0.3*logPval.max()])

    ax.set_xlabel(" -log (pvalue) ")
    ax.set_ylabel(" median XlogP of clusters ")

    plt.show()

    '''
    ax = sns.scatterplot(x="xlogp", y="logPval",
                         hue="altratio", size="csize",
                         sizes=(100, 1000), hue_norm=(0, 7),
                         legend=False, data=clustdf_alt_impact)

    '''


def ont_graph(clusterdf_all, obo_file='ChemOnt_2_1/ChemOnt_2_1.obo',
              entity_names='ChemOnt_2_1/entity_names.txt', out_file=''):
    graph = obonet.read_obo(obo_file)
    id_to_name = {id_: data.get('name') for id_, data in graph.nodes(data=True)}
    nms = pd.read_csv(entity_names, sep='\t', header=None)
    clustdf = clusterdf_all[clusterdf_all['pvalues'] != 1 ]
    clustdf_alt_impact = clustdf[(clustdf['pvalues'] < 0.05) &
                                  (clustdf['csize'] > 1) &
                                  (clustdf['alteredMetabolites'] > 1)]
    clustdf_alt_impact = pd.merge(clustdf_alt_impact, nms, left_on='name', right_on=1)

    clustdf_alt_impact.sort_values(['xlogp'], inplace=True)
    clustdf_alt_impact['order'] = range(clustdf_alt_impact.shape[0])
    clustdf_alt_impact['logPval'] = [-np.log(x) for x in clustdf_alt_impact['pvalues']]

    cmap = {obsv_id:clustdf_alt_impact.loc[obsv_id, "upratio"] for obsv_id in clustdf_alt_impact.index}

    clist = list(sns.color_palette("RdBu_r", 7).as_hex())
    #clist.reverse()
    sm = ScalarMappable(norm=Normalize(vmin=min(list(cmap.values())),
                                       vmax=max(list(cmap.values()))),
                        cmap=ListedColormap(clist))

    s = [20*n for n in clustdf_alt_impact["logPval"]]
    c = [sm.to_rgba(cmap[obsv_id]) for obsv_id in
           clustdf_alt_impact.index]

    labels = {}
    ex = []
    for i in clustdf_alt_impact.index:
        if clustdf_alt_impact.loc[i,0] in id_to_name.keys():
            labels[clustdf_alt_impact.loc[i,0]] = id_to_name[clustdf_alt_impact.loc[i,0]]
        else:
           ex.append(i)

    altpos = {}
    for i in clustdf_alt_impact.index:
        altpos[clustdf_alt_impact.loc[i,0]] = np.array([clustdf_alt_impact.loc[i, "xlogp"], clustdf_alt_impact.loc[i, "logPval"]])

    subgraph = graph.subgraph(labels.keys())
    #pos=nx.fruchterman_reingold_layout(subgraph)
    #pos=nx.circular_layout(subgraph)
    #pos=nx.random_layout(subgraph)
    #pos=nx.spectral_layout(subgraph)
    pos=nx.spring_layout(subgraph)
    #pos = altpos
    nx.draw_networkx_nodes(subgraph,pos,
                           node_color=c,
                           node_size=s,
                           alpha=0.8)
    nx.draw_networkx_edges(subgraph,pos,width=1.0,alpha=0.5)
    nx.draw_networkx_labels(subgraph,pos,labels,font_size=6)
    plt.show()
    if out_file!='':
        G=nx.Graph()
        G.add_nodes_from(subgraph.node)
        G.add_edges_from(list(subgraph.edges.data()))
        nx.set_node_attributes(G, labels, 'labels')
        c = {clustdf_alt_impact.loc[obsv_id,0]:mcolors.to_hex(sm.to_rgba(cmap[obsv_id]))
             for obsv_id in clustdf_alt_impact.index}
        nx.set_node_attributes(G, c, 'color')
        s = {clustdf_alt_impact.loc[i,0]:20*clustdf_alt_impact.loc[i,"logPval"]
             for i in clustdf_alt_impact.index}
        nx.set_node_attributes(G, s, 'size')
        nx.write_gml(G, out_file)



def db2set(entity_names, class_set, db,
           field='InChI', classyfire_name='class_name',
           outfile=None):
    """ Creates class dictionaries from databases containing
        structural information
    Parameters
    ----------
    entity_names : str or pandas.DataFrame
       ClassyFire class to id relationship
    class_set : str or dict
       Reference class dictionary
    db : str or pandas.DataFrame
       Database containing InChI or SMILES and ClassyFire name
    field : str
       InChI or SMILES
    classyfire_name : str
       ClassyFire name present in the database
    outfile : str
       Output file for new set
    Returns
    -------
    class_dict : dict
       Dictionary with keys from ClassyFire ids and
       values with from InChIKeys
    .json : file
       gseapy object
    """
    if type(entity_names)==str:
        entity_names = pd.read_csv(entity_names, sep='\t', header=None)

    if type(class_set)==str:
        with open(class_set) as f:
            class_set = json.load(f)

    if type(db)==str:
        db = pd.read_csv(db, sep='\t')

    if field=='InChI':
        db['InChIKey'] = [Chem.InchiToInchiKey(x)
                          if type(x)==str else '' for x in db['InChI']]
    elif field=='SMILES':
        db['InChIKey'] = [Chem.MolToInchiKey(Chem.MolFromSmiles(x))
                          if type(x)==str else '' for x in db['SMILES']]

    db = pd.merge(db, entity_names, left_on=classyfire_name, right_on=1)

    db2set = {}
    for k in list(class_set.keys()):
        set_list = db.loc[db[0]==k, 'InChIKey'].unique().tolist()
        if len(set_list)>2:
            db2set[k] = set_list

    if outfile!=None:
        with open(outfile, 'w+') as f:
            json.dump(db2set, f, indent=4)

    return db2set
