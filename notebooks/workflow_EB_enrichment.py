"""Enrichment analysis for the GNPS2 EB workflow.

Mirrors notebooks/workflow_simplified.ipynb for the EB task
049e95cf415d4a8fa1a2437635d062d2, comparing the ATTRIBUTE_dx_full metadata
field (POS against NEG).

Run it inside the pyclassrich container:

    docker run --rm -v <repo>:/repo -w /repo -e PYTHONPATH=/repo \
        -e MPLBACKEND=Agg pyclassrich-worker \
        /opt/conda/envs/pyclassrich/bin/python notebooks/workflow_EB_enrichment.py
"""
import os

import matplotlib
matplotlib.use('Agg')

import pandas as pd

from pyclassrich.gnps import Proteosafe
from pyclassrich.models import class_enrichment, impact_plot
from pyclassrich.stats import univariate
from pyclassrich.utils import get_db

TASKID = '049e95cf415d4a8fa1a2437635d062d2'
OUTDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'eb_enrichment')

params = {
    "mzmine": False,
    "mzmine_batch": False,
    "desc": 'description',
    "gnps_taskid": TASKID,
    "gnps_workflow": 'EB',
    "gnps_nap": False,
    "nap_taskid": '',
    "classify": True,
    "normalization": {
        "perform": True,
        "type": 'TIC'
    },
    "enrich": True,
    "comparison": {
        "classes": 'POS,NEG',
        "field": 'ATTRIBUTE_dx_full',
        "test": 'ttest',
        "vthr": 0.5
    },
    "summarize": True,
    "pcoa": True,
    "pcoa_metric": 'canberra',
    "pcoa_norm": True,
    "pcoa_scale": True,
    "email": 'email@gmail.com',
    "chw": '',
    "type": 'canopus',
    "canopus_file": ''
}


def load_compound_classes():
    """ ClassyFire ranks used to name the chemical classes of the EB hits.
    The COCONUT table ships with the pyclassrich container, the sqlite
    compound database is used when the table is not available.
    Returns
    -------
    classes : pandas.DataFrame
        InChI keyed class annotation table.
    """
    psv = os.environ.get('COCONUT_PSV', '/formatdb_flask/COCONUT.psv')
    usecols = ['InChI', 'Identifier', 'kingdom_name', 'superclass_name',
               'class_name', 'subclass_name']
    if os.path.exists(psv):
        db = pd.read_csv(psv, sep='|', usecols=usecols)
    else:
        db = get_db('HMDB,GNPS,DNP,CHEBI,SUPNAT,DRUGBANK,MarinLit,FooDB,'
                    'NPAtlas,COCONUT')
        db = db[usecols]
    return db[~db['InChI'].duplicated()].reset_index(drop=True)


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    # ------------------------------------------------------------------
    # data retrieval
    # ------------------------------------------------------------------
    gnps = Proteosafe(TASKID, 'EB')
    gnps.get_gnps()

    feat = gnps.feat
    meta = gnps.meta
    meta.filename = meta.filename.str.replace(' Peak area', '')

    sampleid = [c for c in feat.columns if ' Peak area' in c]
    samplekey = pd.Series(sampleid).str.replace(r'(.+)\.mzX?ML .+', r'\1', regex=True)
    metakey = meta['filename'].str.replace(r'\.mzXML|\.mzML', '', regex=True)

    print('== EB data retrieval (task %s) ==' % TASKID)
    print('nodes       : %6d rows x %3d columns [%s]'
          % (gnps.gnps.shape + (', '.join(gnps.gnps.columns),)))
    print('annotations : %6d rows, %d #Scan# values outside the feature table'
          % (gnps.dbmatch.shape[0],
             (~gnps.dbmatch['#Scan#'].isin(gnps.gnps['cluster index'])).sum()))
    print('edges       : %6d rows, %d components'
          % (gnps.net.shape[0], gnps.net['ComponentIndex'].nunique()))
    print('features    : %6d rows x %3d columns, %d quantified samples'
          % (feat.shape[0], feat.shape[1], len(sampleid)))
    print('metadata    : %6d rows x %3d columns, %d samples matched'
          % (meta.shape[0], meta.shape[1], len(set(samplekey) & set(metakey))))
    print('spectra     : %6d entries' % len(gnps.spectra))
    print('classes     : %s' % meta['ATTRIBUTE_dx_full'].value_counts().to_dict())

    assert meta['ATTRIBUTE_dx_full'].isin(
        params['comparison']['classes'].split(',')).sum() > 0
    assert set(samplekey) == set(metakey), 'feature and metadata sample names differ'
    assert gnps.dbmatch['#Scan#'].isin(gnps.gnps['cluster index']).all()

    # ------------------------------------------------------------------
    # annotation table, as built in workflow_simplified.ipynb
    # ------------------------------------------------------------------
    tabgnps = pd.merge(gnps.gnps, gnps.dbmatch, left_on='cluster index',
                       right_on='#Scan#', how='left')
    annotated = tabgnps[['cluster index', 'SpectrumID', 'Smiles', 'INCHI',
                         'InChIKey', 'superclass', 'class', 'subclass']]
    annotated = annotated.rename(columns={'cluster index': 'cluster.index',
                                          'superclass': 'superclass_name',
                                          'class': 'class_name',
                                          'subclass': 'subclass_name',
                                          'InChIKey': 'InchiKey'})

    annotated['Identifier'] = annotated['InChI'] = annotated['superclass_id'] = ''
    annotated['class_id'] = annotated['subclass_id'] = annotated['conflict'] = ''
    annotated['ginchikey'] = annotated['kingdom_name'] = ''
    annotated['Score'] = 0

    annotated = annotated[['cluster.index', 'SpectrumID', 'Smiles', 'INCHI',
                           'Identifier', 'Score', 'InChI', 'kingdom_name',
                           'superclass_name', 'class_name', 'subclass_name',
                           'superclass_id', 'class_id', 'subclass_id',
                           'conflict', 'ginchikey', 'InchiKey']]

    tabgnps = tabgnps.rename(columns={'cluster index': 'cluster.index',
                                      'parent mass': 'parent.mass',
                                      'SpectrumID': 'LibraryID'})

    print('node/annotation merge : %d rows x %d columns, %d library hits'
          % (tabgnps.shape[0], tabgnps.shape[1],
             tabgnps['LibraryID'].notna().sum()))

    # ------------------------------------------------------------------
    # univariate test, POS against NEG
    # ------------------------------------------------------------------
    df = univariate(annotated, tabgnps, params, feat, meta,
                    input_dir=OUTDIR, filename='EB_univariate.tsv')

    print('univariate  : %6d rows, %d features with p < 0.05, %d with FDR < 0.05'
          % (df.shape[0], (df['pval'] < 0.05).sum(), (df['cpval'] < 0.05).sum()))

    # ------------------------------------------------------------------
    # chemical classes (ClassyFire ranks from COCONUT) and enrichment
    # ------------------------------------------------------------------
    classes = load_compound_classes()

    # univariate replaces the raw INCHI of the library hits by the empty
    # classyfire placeholders of the annotated table, so the structure and
    # the library name are taken back from dbmatch
    placeholders = ['Identifier', 'kingdom_name', 'superclass_name',
                    'superclass_id', 'class_name', 'class_id',
                    'subclass_name', 'subclass_id', 'conflict', 'ginchikey',
                    'InChI', 'Score']
    uni = df.drop(columns=[c for c in placeholders if c in df.columns])
    lib = gnps.dbmatch[['#Scan#', 'INCHI', 'Compound_Name']].rename(
        columns={'#Scan#': 'cluster.index', 'INCHI': 'InChI',
                 'Compound_Name': 'library_name'})
    uni = pd.merge(uni, lib, on='cluster.index', how='left')
    uni = pd.merge(uni, classes, on='InChI', how='left')

    uni['class_name'] = uni['class_name'].fillna('')
    uni['Identifier'] = uni['Identifier'].fillna('')
    uni['Compound_Name'] = uni['Identifier'].where(uni['Identifier'] != '',
                                                   uni['library_name'])
    uni.to_csv(os.path.join(OUTDIR, 'EB_uni.tsv'), sep='\t', index=False)

    hits = uni['library_name'].notna()
    chemrich = uni.loc[hits, ['Compound_Name', 'class_name', 'InChI', 'pval',
                              'fchange']]
    chemrich.columns = ['Compound_Name', 'Class', 'InChI', 'pvalue', 'foldchange']
    chemrich.to_csv(os.path.join(OUTDIR, 'EB_chemrich.tsv'), sep='\t', index=False)

    print('classes     : %d library hits, %d with an InChI, %d with a '
          'ClassyFire class (%d distinct classes)'
          % (hits.sum(), chemrich['InChI'].fillna('').str.len().gt(0).sum(),
             chemrich['Class'].astype(bool).sum(),
             chemrich.loc[chemrich['Class'] != '', 'Class'].nunique()))

    clusterdf = class_enrichment(chemrich, cfield='Class', nfield='Compound_Name')
    clusterdf.to_csv(os.path.join(OUTDIR, 'EB_clusterdf.tsv'), sep='\t', index=False)

    print('enrichment  : %d classes with more than two members, %d with p < 0.05'
          % (clusterdf.shape[0], (clusterdf['pvalues'] < 0.05).sum()))
    print(clusterdf.sort_values('pvalues')[
        ['name', 'pvalues', 'adjustedpvalue', 'csize', 'upcount', 'downcount',
         'xlogp']].head(10).to_string(index=False))

    try:
        impact_plot(clusterdf, chemrich, cfield='Class', nfield='Compound_Name',
                    scaling_factor=20, top=10, xshift=0.5, outdir=OUTDIR,
                    filename='EB_impact.pdf',
                    fig_title='EB ATTRIBUTE_dx_full POS vs NEG')
        print('impact plot : %s' % os.path.join(OUTDIR, 'EB_impact.pdf'))
    except Exception as error:
        print('impact plot : skipped (%s)' % error)

    print('results in  : %s' % OUTDIR)
    return clusterdf


if __name__ == '__main__':
    main()
