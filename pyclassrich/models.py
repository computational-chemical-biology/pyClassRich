
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
from rdkit.Chem import Descriptors

from matplotlib.backends.backend_pdf import PdfPages

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

from pyclassrich.gnps import Proteosafe

def apply_gsea(nap_taskid, mol_set, params, annotated,
               feat=None, meta=None, outdir='',
               method='ratio_of_classes', clname='',
               cfield=''):
    """ Applies GSEA to contrast defined groups on FBMN feature tables
    Parameters
    ----------
    nap_taskid : str
       NAP task id
    mol_set : str or dict
       Molecular class dictionary
    outdir : str
       Output directory
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

    if type(feat)==type(None) or \
       type(meta)==type(None):
        gnps = Proteosafe(params['gnps_taskid'], 'MZmine')
        gnps.get_gnps()
        feat = gnps.feat
        meta = gnps.meta

    feat = pd.merge(annotated[['cluster.index', 'InchiKey']], feat, left_on='cluster.index', right_on='row ID')
    feat.drop(['cluster.index', 'row m/z', 'row retention time'], axis=1, inplace=True)
    # check if all inchikeys are calculated
    feat.rename(columns={'InchiKey':'NAME', 'row ID':'DESCRIPTION'}, inplace=True)
    feat.columns = feat.columns.str.replace('(.+)\\.mzX?ML .+', '\\1')
    meta['filename'] = meta['filename'].str.replace("\\.mzXML|\\.mzML", "")
    meta['filename'] = meta['filename'].str.replace('\s+$', '')

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
    feat = feat[~((feat['NAME']=='') | (feat['NAME'].isnull()) |
                 (feat['NAME'])==0)]
    feat['DESCRIPTION'] = 0

    ex = feat[feat.columns[2:]].apply(lambda a: sum(a)==0, axis=1)
    if ex.sum():
        feat = feat[~ex]

    # have its own norm param?
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
            if clname not in c:
                continue
            gsea_dict[c] = {}
            for key,value  in mol_set.items():
                if cfield not in key:
                    continue
                cls_vec_tmp = np.array(cls_vec)
                cls_vec_tmp[cls_vec_tmp!=c] = 'others'
                gs_res = gp.gsea(data=feat, # or data='./P53_resampling_data.txt'
                                 gene_sets=value, # enrichr library names
                                 cls=cls_vec_tmp.tolist(), # './data/P53.cls'
                                 # set permutation_type to phenotype if samples >=15
                                 permutation_type='phenotype',
                                 permutation_num=1000, # reduce number to speed up test
                                 outdir=None,  # do not write output to disk
                                 no_plot=True, # Skip plotting
                                 method=method,
                                 processes=4,
                                 format='png')
                gsea_dict[c][key] = gs_res
        return gsea_dict
    else:
        gsea_dict = {}
        for key,value  in mol_set.items():
            if cfield not in key:
                continue
            gsea_dict[key] = {}
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
            gsea_dict[key] = gs_res
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

def class_enrichment(chemrich, cfield='ChemRICHClusters',
                     nfield='CompoundName', method='kstest',
                     pfield='pvalue', params=None, logp=True):
    """ Applies class enrichment
    Parameters
    ----------
    chemrich : pd.DataFrame
       Dataframe containing identifiers, classes and statistics
    cfield : str
       Column label for chemical class
    nfield : str
       Compound identifiers
    method : str
       Enrichment method, either 'kstest' or 'gsea'
    pfield : str
       p-value field
    params : dict
       If method is gsea, required parameters
    logp : bool
       Whether to calculate logp
    Returns
    -------
    chemrich : pd.DataFrame
       Dataframe with enrichment results
    """

    chemdict = Counter(chemrich[cfield])
    clusterids = np.array(list(chemdict.keys()))
    clusterids = clusterids[np.array(list(chemdict.values()))>2]
    clusterids = clusterids[(clusterids!='nan') & (clusterids!='')]

    if method=='ktest':
        clspvals = [kstest(chemrich, x, cfield=cfield, pfield=pfield) for x in clusterids]
        clspvals = np.array(clspvals)
        # May be not needed
        clspvals[clspvals==0] = 2**(-20)
        clusterdf = pd.DataFrame([clusterids, clspvals], index=['name', 'pvalues']).T
        clusterdf['adjustedpvalue'] = multipletests(clusterdf['pvalues'], method = "fdr_bh")[1]
    elif method=='gsea':
        ndict = Proteosafe(params['nap_taskid'], 'NAP')
        clfy = Chemont('scripts/ChemOnt_2_1/ChemOnt_2_1.obo')
        ndict.get_nap()
        ndb =  ndict.ndb
        ndb = ndb.replace(' ', '')
        dblist = [getdb(x, 'scripts/DB/') for x in ndb.split(',')]
        dblist = [pd.read_csv(x, sep='\t') for x in dblist if x is not None]
        db = pd.concat(dblist, ignore_index=True)
        mol_set = db2set(db,
                        'scripts/ChemOnt_2_1/ChemOnt_2_1.obo')
        gsres = apply_gsea(params['nap_taskid'], mol_set,
                          params, ctabsel,
                          clname=params['clname'].split('_')[0],
                          cfield=cfield)
        cls = params['comparison']['classes']
        if cls=='all' or len(cls.split(','))>2:
            res2d = gsres[params['clname'].split('_')[0]]['%s_name' % cfield].res2d
        else:
            res2d = gsres['%x_name' % cfield].res2d
        # May be not needed
        clusterdf = res2d[['pval', 'fdr']].reset_index()
        clusterdf.columns = ['name', 'pvalues', 'adjustedpvalue']
        clusterdf['name'] = [clfy.id_to_name[x] for x in clusterdf.name]
        clusterdf.loc[clusterdf['pvalues']==0, 'pvalues'] = 2**(-20)
        clusterdf.loc[clusterdf['adjustedpvalue']==0, 'adjustedpvalue'] = 2**(-20)
    else:
        raise ValueError("No such method as %s" % method)

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

    if logp:
        mols = [Chem.MolFromInchi(x) for x in chemrich['InChI']]
        # Flag 0?
        chemrich['xlogp'] = [Descriptors.MolLogP(m) if m!=None else 0                                                                                                                             for m in mols]
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


def impact_plot(clusterdf, chemrich, cfield='ChemRICHClusters',
                nfield='CompoundName', scaling_factor=20,
                outdir='', filename='', fig_title=''):
    """ Impact plot
    Parameters
    ----------
    chemrich : pd.DataFrame
       Dataframe containing identifiers, classes and statistics
    clusterdf : pd.DataFrame
       Dataframe with enrichment results
    cfield : str
       Column label for chemical class
    nfield : str
       Compound identifiers
    scaling_factor : int
       Compound identifiers
    outdir : str
       Output directory
    filename : str
       file name
    fig_title : str
       Figure title
    Returns
    -------
    impact_plot : pdf
    """

    chemrich[nfield] = chemrich[nfield].astype(str)
    Compounds = []
    for x in clusterdf['name']:
        cpds = chemrich.loc[chemrich[cfield]==x, nfield]
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


    #s = [scaling_factor*n for n in clustdf_alt_impact["logPval"]]
    s = [scaling_factor*n for n in clustdf_alt_impact["csize"]]

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

    ax.set_xlabel(" median XlogP of clusters ")
    ax.set_ylabel(" -log (pvalue) ")
    ax.set_title(fig_title)

    if outdir!='':
        with PdfPages(os.path.join(outdir, filename)) as pdf:
            pdf.savefig()
            plt.close()
    else:
        plt.show()

    '''
    ax = sns.scatterplot(x="xlogp", y="logPval",
                         hue="altratio", size="csize",
                         sizes=(100, 1000), hue_norm=(0, 7),
                         legend=False, data=clustdf_alt_impact)

    '''

def get_impact_plot(ctabsel, params, cfield, nfield, method,
                    scaling_factor, hierarchy, outdir,
                    top=0):
    """ Calculates chemical enrichment and generates plots
    Parameters
    ----------
    ctabsel : pd.DataFrame
       Dataframe containing identifiers, classes and statistics
    params : dict
       Dict containing the classes being compared
    cfield : str
       Column label for chemical class
    nfield : str
       Compound identifiers
    method : str
       Enrichment method, either 'kstest' or 'gsea'
    scaling_factor : int
       Compound identifiers
    hierarchy : str
       ClassyFire chemical hierarchy
    outdir : str
       Output directory
    top : int
       Number of top lower p-values to plot
    Returns
    -------
    impact_plot : pdf
    chemrich : tsv
    """
    outdir = os.path.join(outdir, 'chemrich')
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    ctabsel.fillna(0, inplace=True)
    cls = params['comparison']['classes']
    if cls=='all' or len(cls.split(','))>2:
        clnames = ctabsel.columns[ctabsel.columns.str.contains('_foldchange')]
        for clname in clnames:
            ctabsel['fchange'] = ctabsel[clname]
            class_table = ctabsel.loc[ctabsel['fchange']!=0, ["Identifier", "pval", "fchange", hierarchy, "InChI"]]
            class_table = class_table[class_table[hierarchy]!=0]
            class_table.sort_values([hierarchy], inplace=True)
            class_table[0] = list(range(class_table.shape[0]))
            # Column names are hard coded to c and nfield pars
            class_table.columns = ['compound names', 'pvalue', 'foldchange', 'class', 'InChI', 'order']
            class_table.reset_index(drop=True, inplace=True)

            prefix = clname.split('_')[0]
            cfilename='impact_plot_%s_%s.tsv' % (prefix, hierarchy)
            params['clname'] = clname
            chemrich = class_enrichment(class_table, cfield=cfield,
                                        nfield=nfield, method=method,
                                        params=params)
            chemrich.to_csv(os.path.join(outdir, cfilename),
                            sep='\t', index=None)
            if top!=0:
                chemrich = chemrich.sort_values(['pvalues']).iloc[:top,:]
            fig_title = '%s vs Average' % prefix
            filename='impact_plot_%s_%s_%s.pdf' % (method, prefix, hierarchy)
            impact_plot(chemrich, class_table, cfield=cfield,
                        nfield=nfield, scaling_factor=scaling_factor,
                        outdir=outdir, filename=filename,
                        fig_title=fig_title)
    else:
        # Are compound names taken from GNPS?
        class_table = ctabsel.loc[ctabsel['fchange']!=0, ["Identifier", "pval", "fchange", hierarchy, "InChI"]]
        class_table = class_table[class_table[hierarchy]!=0]
        class_table.sort_values([hierarchy], inplace=True)
        class_table[0] = list(range(class_table.shape[0]))
        class_table.columns = ['compound names', 'pvalue', 'foldchange', 'class', 'InChI', 'order']
        class_table.reset_index(drop=True, inplace=True)

        chemrich = class_enrichment(class_table, cfield=cfield,
                                        nfield=nfield, method=method,
                                        params=params)
        chemrich.to_csv(os.path.join(outdir, 'chemrich_%s.tsv' % hierarchy),
                        sep='\t', index=None)
        filename='impact_plot_%s_%s.pdf' % (method, hierarchy)
        impact_plot(chemrich, class_table, cfield=cfield,
                    nfield=nfield, outdir=outdir,
                    scaling_factor=scaling_factor,
                    filename=filename,
                    fig_title=' vs '.join(cls.split(',')))

# To be applied 
def plot_classes_by_group(ctabsel):
    """ Upregulated classes by category
    Parameters
    ----------
    ctabsel : pd.DataFrame
       Dataframe containing identifiers, classes and statistics
    cfield : str
       Column label for chemical class
    outdir : str
       Output directory
    filename : str
       file name
    fig_title : str
       Figure title
    Returns
    -------
    count_plots : pdf
    """
    ctabsel['max_class'] = ctabsel.iloc[:,ctabsel.columns.str.contains('_foldchange')].apply(lambda a: a.idxmax(), axis=1)

    df = pd.DataFrame(Counter(ctabsel.loc[ctabsel['pval']<0.05, 'max_class']), index=[0]).T.reset_index()
    df.columns = ['Class', 'Metabolite Count']
    df['Class'] = df['Class'].str.replace('_foldchange', '')
    ax = df.plot.bar(x='Class', y='Metabolite Count', rot=45)

    with PdfPages(os.path.join(outdir, 'significant_barplot.pdf')) as pdf:
        pdf.savefig()
        plt.close()

    df = ctabsel.groupby([cfield, 'max_class']).size().reset_index(name='counts')
    df = df[df[cfield]!=0]
    mx = ctabsel[cfield].value_counts().sort_values().tail(10).index[:9]
    df = df[df[cfield].isin(mx)]
    df['max_name'] = df['max_name'].str.replace('_foldchange', '')

    ax = sns.barplot(x=cfield, y="counts", hue="max_class", data=df)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
    ax.tick_params(labelsize=5)

    with PdfPages(os.path.join(outdir, 'class_location.pdf')) as pdf:
        pdf.savefig()
        plt.close()



def ont_graph(clusterdf_all, obo_file='scripts/ChemOnt_2_1/ChemOnt_2_1.obo',
              out_file='', show=False):
    clfy = Chemont(obo_file)
    id_to_name = clfy.id_to_name
    nms = pd.DataFrame(id_to_name, index=[0]).T.reset_index()
    nms.columns = [0, 1]
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

    subgraph = clfy.graph.subgraph(labels.keys())
    #pos=nx.fruchterman_reingold_layout(subgraph)
    #pos=nx.circular_layout(subgraph)
    #pos=nx.random_layout(subgraph)
    #pos=nx.spectral_layout(subgraph)
    pos = nx.spring_layout(subgraph)
    #pos = altpos
    if show:
        nx.draw_networkx_nodes(subgraph,pos,
                               node_color=c,
                               node_size=s,
                               alpha=0.8)
        nx.draw_networkx_edges(subgraph,pos,width=1.0,alpha=0.5)
        nx.draw_networkx_labels(subgraph,pos,labels,font_size=6)
        plt.show()
    if out_file!='':
        G = nx.Graph()
        G.add_nodes_from(subgraph.nodes)
        G.add_edges_from(list(subgraph.edges.data()))
        nx.set_node_attributes(G, labels, 'labels')
        c = {clustdf_alt_impact.loc[obsv_id,0]:mcolors.to_hex(sm.to_rgba(cmap[obsv_id]))
             for obsv_id in clustdf_alt_impact.index}
        nx.set_node_attributes(G, c, 'color')
        s = {clustdf_alt_impact.loc[i,0]:20*clustdf_alt_impact.loc[i,"logPval"]
             for i in clustdf_alt_impact.index}
        nx.set_node_attributes(G, s, 'size')
        nx.write_gml(G, out_file)

def get_ont_graph(ctabsel, params, cfield, nfield, method,
                  outdir, scaling_factor=10):
    """ Calculates chemical enrichment and generates plots
    Parameters
    ----------
    ctabsel : pd.DataFrame
       Dataframe containing identifiers, classes and statistics
    params : dict
       Dict containing the classes being compared
    cfield : str
       Column label for chemical class
    nfield : str
       Compound identifiers
    method : str
       Enrichment method, either 'kstest' or 'gsea'
    scaling_factor : int
       Compound identifiers
    outdir : str
       Output directory
    Returns
    -------
    ont_graph : gml
    """
    outdir = os.path.join(outdir, 'chemrich')
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    ctabsel.fillna(0, inplace=True)
    cls = params['comparison']['classes']
    if cls=='all' or len(cls.split(','))>2:
        clnames = ctabsel.columns[ctabsel.columns.str.contains('_foldchange')]
        for clname in clnames:
            ctabsel['fchange'] = ctabsel[clname]
            chemlist = []
            for hierarchy  in ['superclass_name', 'class_name',
            'subclass_name']:
                cfield = hierarchy.split('_')[0]
                class_table = ctabsel.loc[ctabsel['fchange']!=0, ["Identifier", "pval", "fchange", hierarchy, "InChI"]]
                class_table = class_table[class_table[hierarchy]!=0]
                class_table.sort_values([hierarchy], inplace=True)
                class_table[0] = list(range(class_table.shape[0]))
                # Column names are hard coded to c and nfield pars
                class_table.columns = ['compound names', 'pvalue',
                'foldchange', cfield, 'InChI', 'order']
                class_table.reset_index(drop=True, inplace=True)

                prefix = clname.split('_')[0]
                cfilename='impact_plot_%s_%s.tsv' % (prefix, hierarchy)
                params['clname'] = clname
                chemrich = class_enrichment(class_table, cfield=cfield,
                                            nfield=nfield, method=method,
                                            params=params)
                chemlist.append(chemrich)
            clusterdf_all = pd.concat(chemlist).reset_index(drop=True)
            filename = 'ont_graph_%s_%s.gml' % (method, clname)
            filename = os.path.join(outdir, filename)
            ont_graph(clusterdf_all, out_file=filename)
    else:
        chemlist = []
        for hierarchy  in ['superclass_name', 'class_name',
        'subclass_name']:
        # Are compound names taken from GNPS?
            class_table = ctabsel.loc[ctabsel['fchange']!=0, ["Identifier", "pval", "fchange", hierarchy, "InChI"]]
            cfield = hierarchy.split('_')[0]
            class_table = class_table[class_table[hierarchy]!=0]
            class_table.sort_values([hierarchy], inplace=True)
            class_table[0] = list(range(class_table.shape[0]))
            class_table.columns = ['compound names', 'pvalue', 'foldchange',
            cfield, 'InChI', 'order']
            class_table.reset_index(drop=True, inplace=True)

            chemrich = class_enrichment(class_table, cfield=cfield,
                                            nfield=nfield, method=method,
                                            params=params)
            chemlist.append(chemrich)

        filename = 'ont_graph_%s.gml' % method
        filename = os.path.join(outdir, filename)
        ont_graph(clusterdf_all, out_file=filename)




def db2set(db, obo_file, field='InChI',
           outfile=None):
    """ Creates class dictionaries from databases containing
        structural information
    Parameters
    ----------
    db : str or pandas.DataFrame
       Database containing InChI or SMILES and ClassyFire name
    obo_file : str
       ChemOnt file path
    field : str
       InChI or SMILES
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
    clfy = Chemont(obo_file)
    clfy.class_hierarchy()
    # There are a few instances where the id is in a different
    # hierarchy
    #set(ctabsel.loc[ctabsel['class_id']!=0, 'class_id'])-set(clfy.classes)

    entity_names = pd.DataFrame(clfy.id_to_name, index=[0]).T
    entity_names.reset_index(inplace=True)
    entity_names.columns = [0, 1]

    if type(db)==str:
        db = pd.read_csv(db, sep='\t')

    # use first block of inchi?
    if field=='InChI':
        db['InChIKey'] = [Chem.InchiToInchiKey(x)
                          if type(x)==str else '' for x in db['InChI']]
    elif field=='SMILES':
        db['InChIKey'] = [Chem.MolToInchiKey(Chem.MolFromSmiles(x))
                          if type(x)==str else '' for x in db['SMILES']]

    class_set = {}
    class_set['superclass_name'] = clfy.supclass
    class_set['class_name'] = clfy.classes
    class_set['subclass_name'] = clfy.subclass

    db2set = {}
    for key, value in class_set.items():
        db_tmp = pd.merge(db, entity_names, left_on=key, right_on=1)
        db2set[key] = {}
        for k in value:
            set_list = db_tmp.loc[db_tmp[0]==k, 'InChIKey'].unique().tolist()
            if len(set_list)>2:
                db2set[key][k] = set_list

    if outfile!=None:
        with open(outfile, 'w+') as f:
            json.dump(db2set, f, indent=4)

    return db2set

class Chemont:
    def __init__(self, chemont_file):
        """ Builds ChemOnt object
        Returns
        -------
        obonet  : obonet DAG
        """
        self.graph = obonet.read_obo(chemont_file)
        self.id_to_name = {id_: data.get('name') for id_, data in self.graph.nodes(data=True)}

    def is_directed_acyclic_graph(self):
        return nx.is_directed_acyclic_graph(self.graph)

    def descendants(self, nodeid):
        return nx.descendants(self.graph, nodeid)

    def ancestors(self, nodeid):
        return nx.ancestors(self.graph, nodeid)

    def edges_to_parents(self, nodeid):
        # Find edges to parent terms
        for child, parent, key in self.graph.out_edges(nodeid, keys=True):
            print(f'• {self.id_to_name[child]} ⟶ {key} ⟶ {self.id_to_name[parent]}')

    def edges_to_children(self, nodeid):
        # Find edges to children terms
        for parent, child, key in self.graph.in_edges(nodeid, keys=True):
            print(f'• {self.id_to_name[child]} ⟵ {key} ⟵ {self.id_to_name[parent]}')

    def class_hierarchy(self):
        # super class
        # 
        node = 'CHEMONTID:0000000'
        supclass = []
        classes = []
        subclass = []
        for parent, child, key in self.graph.in_edges(node, keys=True):
            supclass.append(parent)
            for pt, cd, ky in self.graph.in_edges(parent, keys=True):
                classes.append(pt)
                for p, c, k in self.graph.in_edges(pt, keys=True):
                    subclass.append(p)
        self.supclass = supclass
        self.classes = classes
        self.subclass = subclass

