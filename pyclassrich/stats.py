import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import xmltodict
import requests
import io
import os
import yaml
import json

import qiime2
from qiime2 import Artifact
import pandas as pd
from qiime2.plugins import metadata, feature_table, diversity, emperor
from scipy.spatial.distance import squareform, pdist
import os
import sys
import yaml
import requests
import xmltodict
import json
import io
import re
import skbio

from pyclassrich.gnps import Proteosafe
from pyclassrich.plots import plot_volcano

def multi_class(data, cn, features, method):
    if len(data.groupby(cn)) <= 2:
        raise Exception('ANOVA requires a secondary index with three or more values')
    if method=='anova':
        func = stats.f_oneway
    elif method=='kruskal':
        func = stats.kruskal
    vstats = []
    for col in fnames:
        cls = []
        for k, v in data[[cn,col]].groupby(cn):
            cls.append(v[col])
        vstats.append(func(*cls))
    return pd.DataFrame(vstats)

def univariate(annotation, tabgnps, params,
               feat, meta, nap_taskid='', volcano=True,
               input_dir='', filename=''):
    """ Performs the selected univariate test in all features.
    Parameters
    ----------
    annotation : str or pd.DataFrame
       Class annotation or file path
    tabgnps : str or pd.DataFrame
       Nodes from GNPS' connected components
    params : str or dict
       Params dict or file path
    feat : str or pd.DataFrame
       Feature quantification or file path
    meta : str or pd.DataFrame
       Feature metadata
    nap_taskid : str
       NAP task id
    volcano : bool
       Whether to plot
    input_dir : str
       Input dir
    filename : str
      Output filename
    Returns
    -------
    .tsv : file
       Annotated GNPS/NAP table with statistics
    .pdf : file
       Volcano plot
    """
    if type(annotation)==str:
        doutput = os.path.join(input_dir, annotation)
        ftabgnps = pd.read_csv(doutput, sep='\t')
    else:
        ftabgnps = annotation
    ctab = ftabgnps.groupby('cluster.index').apply(lambda a: a.reset_index(drop=True).loc[0,:])
    ctab.reset_index(drop=True, inplace=True)
    ctab.sort_values(['cluster.index'], inplace=True)
    ctabsel = ctab[['cluster.index', 'Identifier', 'InchiKey',
                    'kingdom_name', 'superclass_name', 'superclass_id',
                    'class_name','class_id', 'subclass_name',
                    'subclass_id', 'Smiles', 'InChI']]

    if type(feat)==str:
        doutput = os.path.join(input_dir, feat)
        feat = pd.read_csv(doutput)
    elif nap_taskid!='':
        ndict = Proteosafe(nap_taskid, 'NAP')
        ndict.get_gnps_taskid()
        gtaskid = ndict.gtaskid
        quant = Proteosafe(gtaskid, 'MZmine')
        quant.get_gnps()
        feat = quant.feat

    if type(meta)==str:
        doutput = os.path.join(input_dir, meta)
        meta = pd.read_csv(doutput, sep='\t')
    elif nap_taskid!='':
        ndict = Proteosafe(nap_taskid, 'NAP')
        ndict.get_gnps_taskid()
        gtaskid = ndict.gtaskid
        quant = Proteosafe(gtaskid, 'MZmine')
        quant.get_gnps()
        meta = quant.meta

    if type(tabgnps)==str:
        doutput = os.path.join(input_dir, tabgnps)
        tabgnps = pd.read_csv(doutput, sep='\t')
    elif nap_taskid!='':
        ndict = Proteosafe(nap_taskid, 'NAP')
        ndict.get_nap()
        tabgnps = ndict.tabgnps

    ctabsel = pd.merge(tabgnps[['cluster.index','parent.mass', 'RTMean', 'LibraryID']],
                       ctabsel, left_on='cluster.index', right_on='cluster.index')

    if type(params)==str:
        doutput = os.path.join(input_dir, params)
        with open(doutput, 'r') as stream:
            params = yaml.safe_load(stream)

    sampleid = pd.DataFrame(feat.columns[3:].tolist(),
                            columns=['sampleid'])

    fnames = feat['row ID'].tolist()
    feat2 = pd.concat([sampleid,
                       feat[feat.columns[3:]].T.reset_index(drop=True)],
                       axis=1)
    feat2.columns = ['sampleid']+fnames
    feat2 = feat2[feat2['sampleid'].str.contains('Peak area')]
    feat2['sampleid'] = feat2['sampleid'].str.replace('(.+)\\.mzX?ML .+', '\\1')
    mnfeat = feat2[feat2.columns[1:]][feat2[feat2.columns[1:]]!=0].min().min()
    # Make it optional?
    feat2.replace(0, mnfeat*(2/3), inplace=True)
    if params['normalization']['perform'] and \
       params['normalization']['type']=='TIC':
        feat2[feat2.columns[1:]] = feat2[feat2.columns[1:]].apply(lambda a: a/sum(a), axis=1)

    meta['filename'] = meta['filename'].str.replace("\\.mzXML|\\.mzML", "")
    # should transformations be applied?
    #feat2[feat2.columns[1:]] = feat2[feat2.columns[1:]].transform(np.log10)
    meta['filename'] = meta['filename'].str.replace('\s+$', '')
    feat2['sampleid'] = feat2['sampleid'].str.replace('\s+$', '')
    feat2 = pd.merge(meta, feat2, left_on='filename',
                     right_on='sampleid')

    if 'filter' in params.keys():
        if params['filter']['perform']:
            fvalues = params['filter']['value'].split(',')
            meta = meta[~meta[params['filter']['field']].isin(fvalues)]

    if 'keep' in params.keys():
        if params['keep']['perform']:
            fvalues = params['keep']['value'].split(',')
            meta = meta[meta[params['keep']['field']].isin(fvalues)]

    if 'filter' in params.keys() or 'keep' in params.keys():
        feat2 = feat2[feat2['sampleid'].isin(meta['filename'])]

    if params['comparison']['test']=='ttest':
        cs = params['comparison']['classes'].split(',')
        cn = params['comparison']['field']
        idx0 = feat2.index[np.where(feat2[cn]==cs[0])]
        idx1 = feat2.index[np.where(feat2[cn]==cs[1])]
        ttest = feat2[fnames].apply(lambda a: stats.ttest_ind(a[idx0], a[idx1]))
        ttest = pd.DataFrame(ttest.tolist())
        pvals = ttest['pvalue']
        pcor = multipletests(pvals, method='fdr_bh')[1]
    elif params['comparison']['test']=='wilcox':
        cs = params['comparison']['classes'].split(',')
        cn = params['comparison']['field']
        idx0 = feat2.index[np.where(feat2[cn]==cs[0])]
        idx1 = feat2.index[np.where(feat2[cn]==cs[1])]
        tkruskal = feat2[fnames].apply(lambda a: stats.kruskal(a[idx0], a[idx1]))
        tkruskal = pd.DataFrame(tkruskal.tolist())
        pvals = tkruskal['pvalue']
        pcor = multipletests(pvals, method='fdr_bh')[1]
    elif params['comparison']['test']=='anova':
        cn = params['comparison']['field']
        if params['comparison']['classes']=='all':
            cs = feat2[cn].value_counts()
            cs = cs[cs>2].index.tolist()
        else:
            cs = params['comparison']['classes'].split(',')
        feat3 = feat2[feat2[cn].isin(cs)]
        tanova = multi_class(feat3, cn, fnames, 'anova')
        pvals = tanova['pvalue']
        pcor = multipletests(pvals, method='fdr_bh')[1]
    elif params['comparison']['test']=='kruskal':
        cn = params['comparison']['field']
        if params['comparison']['classes']=='all':
            cs = feat2[cn].value_counts()
            cs = cs[cs>2].index.tolist()
        else:
            cs = params['comparison']['classes'].split(',')
        feat3 = feat2[feat2[cn].isin(cs)]
        tkruskal = multi_class(feat3, cn, fnames, 'kruskal')
        pvals = tkruskal['pvalue']
        pcor = multipletests(pvals, method='fdr_bh')[1]

    ctabsel = pd.merge(ctabsel,
                       feat[['row ID',  'row m/z',
                             'row retention time']],
                       left_on='cluster.index',
                       right_on='row ID', how='outer')

    if len(cs)==2:
        idx0 = feat2.index[np.where(feat2[cn]==cs[0])]
        idx1 = feat2.index[np.where(feat2[cn]==cs[1])]
        fchange = feat2[fnames].apply(lambda a:
                                      np.mean(a[idx0])/np.mean(a[idx1]))
        if volcano:
            vdf = pd.DataFrame([np.log2(fchange.tolist()),
                                -np.log10(pvals.tolist()),
                                -np.log10(pcor.tolist())]).T
            vdf.columns = ['fold_change', 'p-value', 'fdr']
            vthr = params['comparison']['vthr']
            pthr = -np.log10(0.05)
            vdf['Class'] = ''
            vdf['Class'] = vdf.apply(lambda a:
                                    cs[0] if a['fold_change']>vthr and a['p-value']>pthr
                                    else (cs[1]
                                    if a['fold_change'] < -vthr and a['p-value']>pthr
                                    else 'Below Threshold') , axis=1)
            palette = {cs[0]: 'red',
                       cs[1]: 'blue',
                       'Below Threshold': 'grey'}
            plot_volcano(vdf, mode='pdf',
                         filename=os.path.join(input_dir, 'ggplot_volcano.pdf'))
        ctabsel['fchange'] = fchange
        ctabsel['pval'] =  pvals
        ctabsel['cpval'] = pcor
    else:
        fchlist = []
        for fch in range(len(cs)):
            idx0 = feat2.index[np.where(feat2[cn]==cs[fch])]
            fchlist.append(feat2[fnames].apply(lambda a:
                                          np.mean(a[idx0])/np.mean(a))
                          )
        fchange = pd.concat(fchlist, axis=1)
        fchange.columns = [c+'_foldchange' for c in cs]

        mxfchange = fchange.apply(lambda a: a[np.where(abs(np.log2(a))==np.amax(abs(np.log2(a))))[0][0]], axis=1)

        if volcano:
            vdf = pd.DataFrame([np.log2(mxfchange.tolist()),
                                -np.log10(pvals.tolist()),
                                -np.log10(pcor.tolist())]).T
            vdf.columns = ['fold_change', 'p-value', 'fdr']
            vthr = params['comparison']['vthr']
            pthr = -np.log10(0.05)
            vdf['Class'] = ''
            vdf['Class'] = vdf.apply(lambda a:
                                    'Above Average' if a['fold_change']>vthr and a['p-value']>pthr
                                    else ('Below Average'
                                    if a['fold_change'] < -vthr and a['p-value']>pthr
                                    else 'Below Threshold') , axis=1)
            palette = {'Above Average': 'red',
                       'Below Average': 'blue',
                       'Below Threshold': 'grey'}
            plot_volcano(vdf, palette, mode='pdf',
                         filename=os.path.join(input_dir, 'ggplot_volcano.pdf'))
        ctabsel['fchange'] = mxfchange
        ctabsel['pval'] =  pvals
        ctabsel['cpval'] = pcor
        # Check here
        ctabsel = pd.concat([ctabsel, fchange], axis=1)

    if filename!='':
        doutput = os.path.join(input_dir, filename)
        ctabsel.to_csv(os.path.join(doutput),
                       sep='\t', index=None)
    return ctabsel

nonphyl_metrics = {'braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine',
                       'dice', 'euclidean', 'hamming', 'jaccard', 'kulsinski', 'mahalanobis', 'matching',
                       'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath',
                       'sqeuclidean', 'wminkowski', 'yule'}

def plot_emperor(taskid, nap_taskid, metric):
    """ Creates emperor plot to MZmine quantification table.
    Parameters
    ----------
    taskid : str
       Carcara task id
    nap_taskid : str
       NAP task id
    metric : str
      Dissimilarity measure, options available in nonphyl_metrics dict
    Returns
    -------
    .qzv : file
       uncompressed qsv file
    """

    """emperor_plot.py <taskid> <metric> <nap_taskid>"""
    target_url = "http://proteomics2.ucsd.edu/ProteoSAFe/DownloadResultFile?task=%s&block=main&file=params/" % nap_taskid
    r = requests.get(target_url)
    plist = xmltodict.parse(r.text)['parameters']['parameter']
    gtaskid = [x['#text'] for x in plist if x['@name']=='JOBID'][0]

    url_to_metadata = "http://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task=%s&block=main&file=metadata_table/" % gtaskid
    url_to_quantification = "http://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task=%s&block=main&file=quantification_table/" % gtaskid

    dout = os.path.join(OUTDIR, taskid)
    with open(os.path.join(dout, 'params.yml'), 'r') as stream:
        params = yaml.safe_load(stream)

    sample_metadata = pd.read_csv(io.StringIO(requests.get(url_to_metadata).text), sep='\t')
    #sample_metadata = pd.read_csv(os.path.join(dout,'metadata.tsv'), sep='\t')
    sample_metadata.rename(index=str, columns={"filename": "#SampleID"}, inplace=True)
    sample_metadata.columns = sample_metadata.columns.str.replace('\s', '_')

    if 'filter' in params.keys():
        if params['filter']['perform']:
            fvalues = params['filter']['value'].split(',')
            sample_metadata = sample_metadata[~sample_metadata[params['filter']['field']].isin(fvalues)]

    if 'keep' in params.keys():
        if params['keep']['perform']:
            fvalues = params['keep']['value'].split(',')
            sample_metadata = sample_metadata[sample_metadata[params['keep']['field']].isin(fvalues)]

    sample_metadata.index = sample_metadata['#SampleID']
    sample_metadata.drop(['#SampleID'], axis=1, inplace=True)
    qsample_metadata = qiime2.metadata.Metadata(sample_metadata.apply(lambda a: a.str.replace(' $', '')))

    df = pd.read_csv(io.StringIO(requests.get(url_to_quantification).text))
    df2 = df[df.columns[df.columns.str.contains(' Peak area')]]
    df2.columns = [re.sub('(.+\.mzX?ML) .+', '\\1', a) for a in df2.columns]
    df2.index = df['row ID'].astype(str)
    df2 = df2.T
    if 'filter' in params.keys() or 'keep' in params.keys():
        df2 = df2[df2.index.isin(sample_metadata.index)]
        df2 = df2[df2.columns[df2.apply(lambda a: sum(a)!=0).tolist()]]
    feat = Artifact.import_data("FeatureTable[Frequency]", df2)
    if params['pcoa_norm']:
        df2 = df2.apply(lambda a: a/sum(a), axis=1)

    if params['pcoa_scale']:
        df2 = (df2-df2.mean())/df2.std()

    dm1 = squareform(pdist(df2, metric=metric))
    dm1 = skbio.DistanceMatrix(dm1, ids=df2.index.tolist())
    dm1 = Artifact.import_data("DistanceMatrix", dm1)
    pcoa = diversity.methods.pcoa(dm1)
    emperor_plot = emperor.visualizers.plot(pcoa.pcoa,
                                            qsample_metadata)
    emperor_plot.visualization.export_data(os.path.join(dout,'qzv'))


