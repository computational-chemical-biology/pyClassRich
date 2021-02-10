import os, re, shutil, sys, datetime, time
import xmltodict, json, io, bz2
import subprocess, requests, ftplib
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import sqlalchemy
import numpy as np
from rdkit import Chem
from collections import defaultdict

#from scripts.classwf.classyfire import *
from pyclassrich.gnps import Proteosafe

def load_sqlite(db_name):
    table_name = "compound"

    engine = sqlalchemy.create_engine("sqlite:///%s" % db_name, execution_options={"sqlite_raw_colnames": True})
    df = pd.read_sql_table(table_name, engine)
    return df

def getdb(dbName, struct_db):
    switcher = {
        'HMDB' : "%s/hmdb.tsv",
        'GNPS' : "%s/gnps.tsv",
        'DNP' : "%s/dnp.tsv",
        'CHEBI' : "%s/chebi.tsv",
        'SUPNAT' : "%s/supernatural.tsv",
        'DRUGBANK' : "%s/drugdb.tsv",
        'MarinLit' : "%s/marinlit.tsv",
        'FooDB' : "%s/fooddb.tsv",
        'NPAtlas' : "%s/npatlas.tsv"
    }
    if dbName in switcher.keys():
        return switcher[dbName] % struct_db
    else:
        return None

def get_db(ndb):
    dbname = {
        'HMDB' : "hmdb",
        'GNPS' : "gnps",
        'DNP' : "dnp",
        'CHEBI' : "chebi",
        'SUPNAT' : "supernatural",
        'DRUGBANK' : "drugdb",
        'MarinLit' : "marinlit",
        'FooDB' : "fooddb",
        'NPAtlas' : "npatlas"
    }
    local_path = f'{os.path.expanduser("~")}/.local/pyclassrich'
    if not os.path.exists(local_path):
        os.makedirs(local_path)

    local_path = f'{os.path.expanduser("~")}/.local/pyclassrich/compounds.db'
    if not os.path.exists(local_path):
        print('Retrieving databse...')
        ftp = ftplib.FTP('seriema.fcfrp.usp.br') #connect to the ftp server
        ftp.login('ftpuser', 'menosfrota')     #using your credentials here

        filedata = open(local_path, 'wb')
        SourceFilename = 'pyclassrich/compounds.db'
        ftp.retrbinary('RETR '+SourceFilename, filedata.write)
        filedata.close()
        ftp.quit()

    print('Loading database...')
    db = load_sqlite(local_path)
    ndb = ndb.split(',')
    ldb = []
    for n in ndb:
        if n in dbname.keys():
            ldb.append(dbname[n])
        else:
            print(f'Database: {n}, not found')
    db = db[db['ref'].isin(ldb)].reset_index(drop=True)

    return db

def inchi2pth(k, db='kegg'):
    if db=='kegg':
        req = requests.get("http://cts.fiehnlab.ucdavis.edu/rest/convert/InChIKey/KEGG/%s" % k)
        if req.status_code!=200:
           return ''
        else:
            return req.json()[0]['results']
    elif db=='biocyc':
        req = requests.get("http://cts.fiehnlab.ucdavis.edu/rest/convert/InChIKey/BioCyc/%s" % k)
        if req.status_code!=200:
           return ''
        else:
            return req.json()[0]['results']

def classify(proteosafe):
    """ Classifies structural annotations from GNPS/NAP
        in decreasing order of confidence
    Parameters
    ----------
    proteosafe : Proteosafe
       Proteosafe object
    Returns
    -------
    annotated : pd.DataFrame
       fully annotated GNPS/NAP table
    """
    proteosafe.get_nap()
    tabgnps = proteosafe.tabgnps
    net = proteosafe.net
    mlist = proteosafe.mlist
    lid = proteosafe.lid
    fusion = proteosafe.fusion
    consensus = proteosafe.consensus
    ndb = proteosafe.ndb
    ndb = ndb.replace(' ', '')

    # Retrieve candidates from in silico fragmentation
    # following the priority order Fusion > Consensus (fusion2) > MetFrag (score)
    sellist = []
    for i in range(tabgnps.shape[0]):
        try:
            tmp = pd.DataFrame(mlist[i])
            ltmp = pd.DataFrame(lid[i])
        except:
            continue
        if tmp.shape[0]==0 and ltmp.shape[0]==0:
            continue
        if 'Identifier' not in tmp.columns:
            pass
        elif len(set(tmp['Identifier']).intersection(set(ltmp['Identifier'])))!=tmp.shape[0]:
            print('Found mismatch')
            break
        if tmp.shape[0] > 3:
            if 'fusion' in tmp.columns:
                tmp['fusion'] = tmp['fusion'].astype(float)
                tmp.sort_values(['fusion'], ascending=False, inplace=True)
                g = tmp[tmp['fusion']==1]['MCSS'].tolist()[0]
                sel = tmp[tmp['MCSS']==g][['Identifier', 'fusion']]
                sel.loc[:,'cluster.index'] = tabgnps.loc[i, 'cluster.index']
                sel.loc[:,'Ncandidates'] = tmp.shape[0]
                sel.rename(columns={'fusion':'Score'}, inplace=True)
            elif 'fusion2' in tmp.columns:
                tmp['fusion2'] = tmp['fusion2'].astype(float)
                tmp.sort_values(['fusion2'], ascending=False, inplace=True)
                g = tmp[tmp['fusion2']==1]['MCSS'].tolist()[0]
                sel = tmp[tmp['MCSS']==g][['Identifier', 'fusion2']]
                sel.loc[:,'cluster.index'] = tabgnps.loc[i, 'cluster.index']
                sel.loc[:,'Ncandidates'] = tmp.shape[0]
                sel.rename(columns={'fusion2':'Score'}, inplace=True)
            elif 'Score' in tmp.columns:
                if not (tmp['Score']==1).any():
                    continue
                tmp['Score'] = tmp['Score'].astype(float)
                tmp.sort_values(['Score'], ascending=False, inplace=True)
                g = tmp[tmp['Score']==1]['MCSS'].tolist()[0]
                sel = tmp[tmp['MCSS']==g][['Identifier', 'Score']]
                sel.loc[:,'cluster.index'] = tabgnps.loc[i, 'cluster.index']
                sel.loc[:,'Ncandidates'] = tmp.shape[0]
        else:
            sel = ltmp[['Identifier', 'Score']]
            sel.loc[:,'cluster.index'] = tabgnps.loc[i, 'cluster.index']
            sel.loc[:,'Ncandidates'] = ltmp.shape[0]
        sellist.append(sel)

    db = get_db(ndb)

    seltab = pd.concat(sellist)
    seltab = pd.merge(seltab, db[['Identifier', 'InChI', 'kingdom_name',
                                  'superclass_name', 'class_name', 'subclass_name',
                                  'superclass_id', 'class_id', 'subclass_id']],
                      left_on='Identifier', right_on='Identifier')

    ftabgnps = pd.merge(tabgnps[['cluster.index', 'SpectrumID', 'Smiles', 'INCHI']], seltab,
                        left_on='cluster.index', right_on='cluster.index', how='left')

    ftabgnps['conflict'] = ''
    ftabgnps['ginchikey'] = ''
    ftabgnps.fillna('', inplace=True)
    for j in range(ftabgnps.shape[0]):
        smi = ftabgnps.loc[j, 'Smiles']
        if type(smi)!=float and smi!='':
            mol = Chem.MolFromSmiles(smi)
            sinchi =  Chem.MolToInchi(mol)
            skey = Chem.InchiToInchiKey(sinchi)
            ftabgnps.loc[j, 'INCHI'] = sinchi
        else:
            skey = None
        inchi = ftabgnps.loc[j, 'INCHI']
        if type(inchi)!=float:
            sinchi = Chem.InchiToInchiKey(inchi)
        if type(smi)!=float and type(inchi)!=float and sinchi!=None and skey!=None:
            skey1 = skey.split('-')[0]
            sinchi1 = sinchi.split('-')[0]
            if skey1==sinchi1:
                ftabgnps.loc[j, 'ginchikey'] = skey
                ftabgnps.loc[j, 'conflict'] = 'no'
            else:
                ftabgnps.loc[j, 'ginchikey'] = sinchi
                ftabgnps.loc[j, 'conflict'] = 'yes'
        elif type(smi)!=float and skey!=None:
            ftabgnps.loc[j, 'ginchikey'] = skey
        elif type(inchi)!=float and sinchi!=None:
            ftabgnps.loc[j, 'ginchikey'] = sinchi
        else:
            continue

    ftabgnps.loc[ftabgnps['INCHI']!='', 'InChI'] = ftabgnps.loc[ftabgnps['INCHI']!='', 'INCHI']
    ftabgnps['InchiKey'] = [Chem.InchiToInchiKey(x) for x in ftabgnps['InChI']]

    inchi = ftabgnps.loc[(ftabgnps['kingdom_name']=='') & (ftabgnps['InChI']!=''), 'InChI']
    if len(inchi):
        cnames = ['superclass_name', 'class_name', 'subclass_name',
                  'superclass_id', 'class_id', 'subclass_id']
        inchi_replace = pd.merge(pd.DataFrame(inchi), db[~db['InChI'].duplicated()], how='left').fillna('')
        ftabgnps.loc[inchi.index, cnames] = inchi_replace.loc[:, cnames]
    return ftabgnps

def getLibra_id_link(i, tab):
    burl = 'https://gnps.ucsd.edu/ProteoSAFe/result.jsp?task='
    purl = '&view=view_all_annotations_DB#%7B%22main.Compound_Name_input%22%3A%22'
    eurl = '%22%7D'
    gid = tab.loc[0,'GNPSLinkout_Cluster'].split('&')[0].split('=')[1]
    if i!='':
        return '<a href=\"'+ burl+gid+purl+str(i)+eurl+ '\" target=\"_blank\">'+str(i)+'</a>'
    else:
        return ''

def get_cluster_spectra_link(i, tab):
    burl = 'https://gnps.ucsd.edu/ProteoSAFe/result.jsp?task='
    purl = '&view=view_all_clusters_withID&show=true#%7B%22main.cluster%20index_lowerinput%22%3A%22'
    eurl = '%22%2C%22main.cluster%20index_upperinput%22%3A%22'
    eeurl = '%22%7D'
    gid = tab.loc[0,'GNPSLinkout_Cluster'].split('&')[0].split('=')[1]
    if i!='':
        return '<a href=\"'+ burl+gid+purl+str(i)+eurl+ str(i)+eeurl+'\"target=\"_blank\">'+'N/A'+'</a>'
    else:
        return ''

def get_library_id_link(i, tab):
    burl = 'https://gnps.ucsd.edu/ProteoSAFe/result.jsp?task='
    purl = '&view=view_all_annotations_DB#%7B%22main.Compound_Name_input%22%3A%22'
    eurl = '%22%7D'
    gid = tab.loc[0,'GNPSLinkout_Cluster'].split('&')[0].split('=')[1]
    if i!='':
        return '<a href=\"'+ burl+gid+purl+str(i)+eurl+ '\" target=\"_blank\">'+str(i)+'</a>'
    else:
        return ''

def linkpattern(id):
    if bool(re.match('^CCMSLIB', id)):
        url = "http://gnps.ucsd.edu/ProteoSAFe/gnpslibraryspectrum.jsp?SpectrumID="+id
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'
    elif bool(re.match('^HMDB', id)):
        url = "http://www.hmdb.ca/metabolites/"+id
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'
    elif bool(re.match('^SN', id)):
        url = "http://bioinf-applied.charite.de/supernatural_new/index.php?site=search5&np_id="+id
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'
    elif bool(re.match('^CHEBI', id)):
        url = "https://www.ebi.ac.uk/chebi/searchId.do?chebiId="+id
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'
    elif bool(re.match('^CID', id)):
        url = "https://pubchem.ncbi.nlm.nih.gov/compound/"+id
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'
    elif bool(re.match('^Chemical-Structure\\.', id)):
        url = "http://www.chemspider.com/"+id+".html"
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'
    elif bool(re.match('^cs', id)):
    #elif bool(re.match('^\d{3,}', id)):
        # provisory for marinlit
        #id2 = id.zfill(12)
        #url = "http://pubs.rsc.org/marinlit/compound/cs"+id2
        url = "http://pubs.rsc.org/marinlit/compound/"+id
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'
    elif bool(re.match('^FDB', id)):
        url = 'http://foodb.ca/compounds/'+id
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'
    elif bool(re.match('^DB', id)):
        url = 'https://www.drugbank.ca/drugs/'+id
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'
    elif bool(re.match('^BGC\d{7,}', id)):
        url = 'https://mibig.secondarymetabolites.org/repository/'+id+'/index.html#cluster-1'
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'
    elif bool(re.match('^NOR\d{5,}', id)):
        url = 'http://bioinfo.lifl.fr/norine/result.jsp?ID='+id
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'
    elif bool(re.match('[A-Z]{3,}', id)):
        # provisory for marinlit
        url = 'http://dnp.chemnetbase.com/faces/chemical/FullScreenEntry.xhtml?product=DNP&fromExport=falsePOSI'+id+'&id='+id
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'
    else:
        url = "https://pubchem.ncbi.nlm.nih.gov/search/#collection=compounds&query_type=text&query="+id
        return '<a href=\"'+ url + '\" target=\"_blank\">'+id+'</a>'

def html_tab(nap_taskid, out_dir,
             fully_annotated_gnps,
             annotated_gnps_nap):
    """ Classifies structural annotations from GNPS/NAP
        in decreasing order of confidence
    Parameters
    ----------
    nap_taskid : str
       NAP task id
    out_dir : str
      Output directory
    fully_annotated_gnps: str or pandas.DataFrame
      DataFrame containing structural annotation
    annotated_gnps_nap: str or pandas.DataFrame
      DataFrame containing structural annotation and univariate results
    Returns
    -------
    .tsv : file
       fully annotated GNPS/NAP table saved to OUTDIR
    stats : pandas.DataFrame
       Formatted DataFrame
    """
    if type(fully_annotated_gnps)==str:
        df = pd.read_csv(os.path.join(out_dir, fully_annotated_gnps), sep='\t')
    else:
        df = fully_annotated_gnps
    if type(annotated_gnps_nap)==str:
        stat = pd.read_csv(os.path.join(out_dir, annotated_gnps_nap), sep='\t')
    else:
        stat = annotated_gnps_nap

    stat.drop(['InchiKey','kingdom_name', 'superclass_name'], axis=1, inplace=True)

    ndict = Proteosafe(nap_taskid, 'NAP')
    ndict.get_nap()

    gnps = Proteosafe(ndict.gtaskid, 'MZmine')
    gnps.get_gnps()
    gnps = gnps.gnps

    stat.fillna('', inplace=True)
    stat['LibraryID'] = stat['LibraryID'].apply(lambda a: getLibraIDlink(a, gnps) if a!='' else '')
    stat.loc[stat['LibraryID']=='', 'LibraryID'] = 'N/A'
    stat.loc[stat['LibraryID']=='N/A', 'LibraryID'] = stat.loc[stat['LibraryID']=='N/A', 'cluster.index'].apply(lambda a: getClusterSpectralink(a, gnps))

    stat['Identifier'] = stat['Identifier'].apply(lambda a: '#'.join([linkpattern(x) for x in list(a.split(','))]) if a!='' else '')

    stat['parent.mass'] = stat['row m/z']
    stat['RTMean'] = stat['row retention time']
    stat.to_csv(os.path.join(out_dir,'tsv2html.txt'), index=None, sep='\t')
    return stat

