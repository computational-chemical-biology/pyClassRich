import xmltodict
import os
import subprocess
import sys
import shutil
import requests
import datetime
import json
import io
import pandas as pd
import numpy as np
import ftplib
import bz2
from rdkit import Chem
import re
from collections import defaultdict
import time

#from scripts.classwf.classyfire import *
from pyclassrich.gnps import Proteosafe


def run_mzmine(mzmine_path, param_file,
               input_dir, output_dir,
               aux_param_file):
    """ Runs MZmine with especified parameters file.
    Parameters
    ----------
    mzmine_path : str
       mzmine sh script path
    param_file : str
       mzmine xml param file
    input_dir : str
       mzML/mzXML input directory
    output_dir : str
       output directory for .csv and .mgf
    aux_param_file : str
       parameter file enforcing some fields
    Returns
    -------
    .csv : file
       Feature quantification table
    .mgf : file
       Fragmentation spectra file
    """
    with open(param_file, 'rb') as f:
        xml = xmltodict.parse(f)

    with open(aux_param_file, 'rb') as f:
        auxxml = xmltodict.parse(f)

    pth = os.path.abspath(input_dir)
    opth = os.path.abspath(output_dir)

    fls = os.listdir(input_dir)
    fls = [os.path.join(pth, x) for x in fls]

    steps = [xml['batch']['batchstep'][i]['@method'] for i in
             range(len(xml['batch']['batchstep']))]

    step0 = 'rawdataimport.RawDataImportModule'
    sfstep = 'csvexport.CSVExportModule'
    fstep_old = 'gnpsexport.GNPSExportAndSubmitModule'
    # For mzmine 2.53 
    fstep = 'gnpsexport.fbmn.GnpsFbmnExportAndSubmitModule'

    if not any(step0 in x for x in steps):
        xml['batch']['batchstep'].insert(0, auxxml['batch']['batchstep'][0])

    if not any(sfstep in x for x in steps):
        xml['batch']['batchstep'].append(auxxml['batch']['batchstep'][8])

    if not any(fstep in x or fstep_old in x for x in steps):
        xml['batch']['batchstep'].append(auxxml['batch']['batchstep'][9])

    for i in range(len(steps)):
        if 'RawDataImportModule' in xml['batch']['batchstep'][i]['@method']:
            xml['batch']['batchstep'][i]['parameter']['file'] = fls
        if type(xml['batch']['batchstep'][i]['parameter'])==list:
            for j in range(len(xml['batch']['batchstep'][i]['parameter'])):
                if 'Raw data files' in xml['batch']['batchstep'][i]['parameter'][j]['@name']:
                    xml['batch']['batchstep'][i]['parameter'][0]['@type'] = 'BATCH_LAST_FILES'
                if 'Mass list' in xml['batch']['batchstep'][i]['parameter'][j]['@name']:
                    xml['batch']['batchstep'][i]['parameter'][j]['#text'] = 'masses'
        if 'CSVExportModule' in xml['batch']['batchstep'][i]['@method']:
            xml['batch']['batchstep'][i]['parameter'][1]['current_file'] = os.path.join(opth, 'feature_table.csv')
            xml['batch']['batchstep'][i]['parameter'][1]['last_file'] = os.path.join(opth, 'feature_table.csv')
        if 'gnpsexport' in xml['batch']['batchstep'][i]['@method']:
            xml['batch']['batchstep'][i]['parameter'][1]['current_file'] = os.path.join(opth, 'spectra.mgf')
            xml['batch']['batchstep'][i]['parameter'][1]['last_file'] = os.path.join(opth, 'spectra.mgf')
        if 'SiriusExport' in xml['batch']['batchstep'][i]['@method']:
            xml['batch']['batchstep'][i]['parameter'][1]['#text'] = os.path.join(opth, 'ms1_ms2.mgf')

    oxml = os.path.join(opth, 'result.xml')
    with open(oxml, 'w+') as result_file:
        result_file.write(xmltodict.unparse(xml, pretty=True))

    subprocess.call([mzmine_path, oxml])

#Base url is just something like massive.ucsd.edu
#parameters is a map
def invoke_workflow(base_url, parameters, login, password):
    username = login
    password = password
    s = requests.Session()
    payload = {
        'user' : username,
        'password' : password,
        'login' : 'Sign in'
    }
    r = s.post('https://' + base_url + '/ProteoSAFe/user/login.jsp', data=payload, verify=False)
    r = s.post('https://' + base_url + '/ProteoSAFe/InvokeTools', data=parameters, verify=False)
    task_id = r.text
    print(task_id)
    if len(task_id) > 4 and len(task_id) < 60:
        print("Launched Task: : " + r.text)
        return task_id
    else:
        return None


def submit_proteosafe(taskid, input_dir, output_dir, invokeParametersFB,
                      invokeParametersNAP, username, password,
                      folderName, email, gnps_taskid=None):
    """ Submits a ProteoSAFe job.
    Parameters
    ----------
    taskid : str
       Carcara task id, or other arbitrary string
    input_dir : str
       Input directory containing spectra, metadata
       and feature table
    output_dir : str
       Output directory
    invokeParametersFB : str
       GNPS FBMN param dict
    invokeParametersNAP : str
       NAP FBMN param dict
    gnps_taskid : str
       parameter file enforcing some fields
    username : str
       GNPS user name
    password : str
       GNPS user password
    folderName : str
       Massive project name
    Returns
    -------
    params.json : file
       parameter file saved to OUTDIR
    paramlist : list
       list of dictionary parameters
    """

    doutput = input_dir
    username = username
    server = 'massive.ucsd.edu'
    #taskid = args[0]

    paramlist = {}
    paramlist['gnps'] = invokeParametersFB
    paramlist['nap'] = invokeParametersNAP

    if gnps_taskid==None:
        ftp_connection = ftplib.FTP(server, username, password, timeout=30)

        folderName = folderName
        if folderName not in ftp_connection.nlst():
        #do needed task 
            ftp_connection.mkd(folderName)

        file = open(os.path.join(doutput, 'spectra.mgf'),'rb')
        ftp_connection.storbinary('STOR '+'/{0}/{1}_spectra.mgf'.format(*[folderName, taskid]), file)
        file.close()                                    # close file and FTP

        file = open(os.path.join(doutput, 'metadata.tsv'),'rb')
        ftp_connection.storbinary('STOR '+'/{0}/{1}_metadata.tsv'.format(*[folderName, taskid]), file)
        file.close()                                    # close file and FTP

        file = open(os.path.join(doutput, 'feature_table.csv'),'rb')
        ftp_connection.storbinary('STOR '+'/{0}/{1}_feature_table.csv'.format(*[folderName, taskid]), file)
        file.close()                                    # close file and FTP

        ftp_connection.quit()

        invokeParametersFB['email'] = email
        invokeParametersFB['metadata_table'] = 'f.{0}/{1}/{2}_metadata.tsv;'.format(*[username, folderName, taskid])
        invokeParametersFB['quantification_table'] = 'f.{0}/{1}/{2}_feature_table.csv;'.format(*[username, folderName, taskid])
        invokeParametersFB['spec_on_server'] = 'f.{0}/{1}/{2}_spectra.mgf;'.format(*[username, folderName, taskid])

        gnps_taskid = invoke_workflow("gnps.ucsd.edu", invokeParametersFB, username, password)

    paramlist['gnps_taskid'] = gnps_taskid
    task_status_url = "https://gnps.ucsd.edu/ProteoSAFe/status_json.jsp?task=%s" % (gnps_taskid)
    task_status = requests.get(task_status_url).json()

    while task_status['status']=='RUNNING':
        task_status = requests.get(task_status_url).json()

    if task_status['status']=='FAILED':
        raise ValueError("GNPS job failed")

    invokeParametersNAP['email'] = email
    invokeParametersNAP['JOBID'] = gnps_taskid

    njobid = invoke_workflow("proteomics2.ucsd.edu", invokeParametersNAP, username, password)
    task_status_url = "https://proteomics2.ucsd.edu/ProteoSAFe/status_json.jsp?task=%s" % (njobid)
    task_status = requests.get(task_status_url).json()

    while task_status['status']=='RUNNING':
        task_status = requests.get(task_status_url).json()

    if task_status['status']=='FAILED':
        raise ValueError("NAP job failed")

    paramlist['nap_taskid'] = njobid

    doutput = output_dir
    with open(os.path.join(doutput, 'params.json'), 'w') as f:
        json.dump(paramlist, f, indent=4, sort_keys=True)

    return paramlist

def getdb(dbName, struct_db):
    switcher = {
        'HMDB' : "%s/hmdb_21_02_2018_formatted.txt",
        'GNPS' : "%s/gnpsCLASS.txt",
        'DNP' : "%s/dnpCLASS.txt",
        'DNP_PUBLIC' : "%s/dnp_public.tsv",
        'CHEBI' : "%s/chebiCLASS.txt",
        'SUPNAT' : "%s/supernaturalCLASS.txt",
        'DRUGBANK' : "%s/drugdb_formatted.txt",
        'MarinLit' : "%s/MarinLit_04_03_2018.txt",
        'FooDB' : "%s/fooddb_04_03_2018.txt",
        'NPAtlas' : "%s/np_atlas_2019_02.tsv"
    }
    if dbName in switcher.keys():
        return switcher[dbName] % struct_db
    else:
        return None

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

def poll(in_process, sleep=10):
    taxon_order = ['kingdom', 'superclass', 'class', 'subclass']
    s = requests.Session()
    records = []
    in_progress_counts = defaultdict(int)
    while in_process:
        id_, query = in_process.pop(0)
        req = s.get("http://classyfire.wishartlab.com/queries/%s.json" % query)
        if req.status_code != 200:
            print("GET something went wrong...")
            continue
        try:
            result = req.json()
        except:
            in_process.append((id_, query))
            in_progress_counts[query] += 1
            continue
        if result['classification_status'] != 'Done':
            if in_progress_counts[query] < 100:
                in_process.append((id_, query))
                in_progress_counts[query] += 1
                continue
        for entry in result['entities']:
            if sum([bool(re.search('report', x)) for x in list(entry.keys())]):
                continue
            if type(entry)==list:
                entry = entry[0]
            if not sum([bool(re.search('identifier', x)) for x in list(entry.keys())]):
                continue
            entry_results = {'#SampleID': entry['identifier']}
            if 'ancestors' in entry:
                for ai, anc in enumerate(entry['ancestors']):
                    anc = anc.replace(' ','_').replace('-','_').replace('/','_').replace(',','_')
                    entry_results[anc] = True
            if 'substituents' in entry:
                for si, subst in enumerate(entry['substituents']):
                    substituent = subst.replace(' ','_').replace('-','_').replace('/','_').replace(',','_')
                    entry_results[substituent] = True
            for taxon in taxon_order:
                if taxon in entry:
                    if entry[taxon] is not None:
                        entry_results['%s_id' % taxon] = entry[taxon].get('chemont_id', None)
                        entry_results['%s_name' % taxon] = entry[taxon].get('name', None)
            #for ci, ch in enumerate(entry['predicted_chebi_terms']):
            #    entry_results['predicted_chebi_term_%d' % ci] = ch
            if 'molecular_framework' in entry:
                entry_results['molecular_framework'] = entry['molecular_framework']
            if 'smiles' in entry:
                entry_results['smiles'] = entry['smiles']
            if 'inchikey' in entry:
                entry_results['inchikey'] = entry['inchikey']
            records.append(entry_results)
        time.sleep(sleep)
    #return pd.DataFrame(records)
    cns = ['inchikey', 'superclass_name', 'superclass_id',
           'class_name', 'class_id', 'subclass_name', 'subclass_id']
    recorddf = pd.DataFrame(records)
    if len(recorddf) and recorddf.columns.isin(cns).sum()==len(cns):
        return recorddf[cns]
    else:
        return None

def get_class(queries, chunksize = 10, sleep=10):
    in_process = []
    s = requests.Session()
    for i in range(0, len(queries), chunksize):
        to_query = queries[i:i+chunksize]
        label = 'id%d' % i
        payload = {'query_input': '\n'.join(to_query), 'query_type': 'STRUCTURE', 'label': label}
        req = s.post('http://classyfire.wishartlab.com/queries.json', data=json.dumps(payload),
                     headers={'Content-Type': "application/json"})
        if req.status_code != 201:
            print("POST Something went wrong...")
            print(req.status_code)
            print(req.content)
        else:
            #in_process.append((label, req.json()['id']))
            in_process.append(poll([(label, req.json()['id'])], sleep))
        time.sleep(sleep)
    return in_process

def classify(output_dir, nap_taskid, struct_db, chemont,
             outfile='', script_path=''):
    """ Classifies structural annotations from GNPS/NAP
        in decreasing order of confidence
    Parameters
    ----------
    output_dir : str
       Output directory
    nap_taskid : str
       NAP task id
    struct_db : str
       Structural DB directory
    chemont : str
       Ontology file from ClassyFire
    outfile : str
       Output file name
    script_path : str
       path to rda2json conversion script, assumes
       Rscript is on the path.
    Returns
    -------
    .tsv : file
       fully annotated GNPS/NAP table saved to output_dir
    status : str
       status code DONE
    """

    nptaskid = nap_taskid
    doutput = output_dir
    # call napviewer or local?
    if script_path=='':
        ndict = Proteosafe(nptaskid, 'NAP')
        ndict.get_nap()
        tabgnps = ndict.tabgnps
        net = ndict.net
        mlist = ndict.mlist
        lid = ndict.lid
        fusion = ndict.fusion
        consensus = ndict.consensus
        ndb =  ndict.ndb
        ndb = ndb.replace(' ', '')
    else:
        subprocess.call(['Rscript', '%s/rda2json2.R' % script_path, nptaskid, doutput])
        tabgnps = pd.read_csv(os.path.join(doutput, 'tabgnps.tsv'), sep='\t')
        net = pd.read_csv(os.path.join(doutput, 'net.tsv'), sep='\t')

        with open(os.path.join(doutput, 'mlist.json'), 'r') as f:
            mlist = json.load(f)

        with open(os.path.join(doutput, 'fusion.json'), 'r') as f:
            fusion = json.load(f)

        with open(os.path.join(doutput, 'lid.json'), 'r') as f:
            lid = json.load(f)

        with open(os.path.join(doutput, 'consensus.json'), 'r') as f:
            consensus = json.load(f)

        with open(os.path.join(doutput, 'params.xml'), 'r') as f:
            nppar = xmltodict.parse(f.read())

        ndb = [x['#text'] for x in nppar['parameters']['parameter'] if x['@name']=='DATABASE'][0]
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

    # check if there is in house database
    dblist = [getdb(x, struct_db) for x in ndb.split(',')]
    dblist = [pd.read_csv(x, sep='\t') for x in dblist if x is not None]
    db = pd.concat(dblist, ignore_index=True)

    seltab = pd.concat(sellist)
    seltab = pd.merge(seltab, db[['Identifier', 'InChI', 'kingdom_name',
                                  'superclass_name', 'class_name', 'subclass_name']],
                      left_on='Identifier', right_on='Identifier')


    ftabgnps = pd.merge(tabgnps[['cluster.index', 'SpectrumID', 'Smiles', 'INCHI']], seltab, left_on='cluster.index', right_on='cluster.index', how='left')

    ftabgnps['conflict'] = ''
    ftabgnps['ginchikey'] = ''
    ftabgnps.fillna('', inplace=True)
    for j in range(ftabgnps.shape[0]):
        smi = ftabgnps.loc[j, 'Smiles']
        if type(smi)!=float and smi!='':
            mol = Chem.MolFromSmiles(smi)
            skey = Chem.MolToInchiKey(mol)
            ftabgnps.loc[j, 'INCHI'] = Chem.MolToInchi(mol)
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

    entity_dict = {}
    entity_names = pd.read_csv(chemont, sep='\t',
                               header=None)

    entity_names.apply(lambda a: entity_dict.update({a[1]:a[0]}), axis=1)
    entity_names.apply(lambda a: entity_dict.update({a[0]:a[1]}), axis=1)

    inchi = ftabgnps.loc[ftabgnps['kingdom_name']=='', 'InChI']
    if len(inchi):
        # improve speed here
        dflist = get_class(inchi.tolist())
        dflist = pd.concat(dflist)
        dflist['inchikey'] = dflist['inchikey'].str.replace('InChIKey=', '')
        dflist.fillna('', inplace=True)
        ftabgnps = pd.merge(ftabgnps,
                          dflist[['inchikey', 'superclass_id','class_id', 'subclass_id']],
                          left_on='InchiKey', right_on='inchikey', how='left')
        ftabgnps.fillna('', inplace=True)
        for idx in ftabgnps.index:
            if ftabgnps.loc[idx, 'superclass_name']=='' and \
               ftabgnps.loc[idx, 'superclass_id']!='':
                ftabgnps.loc[idx, 'superclass_name'] = entity_dict[ftabgnps.loc[idx, 'superclass_id']]
                if ftabgnps.loc[idx,'class_id']!='':
                    ftabgnps.loc[idx, 'class_name'] = entity_dict[ftabgnps.loc[idx,'class_id']]
                if ftabgnps.loc[idx,'subclass_id']!='':
                    ftabgnps.loc[idx, 'subclass_name'] = entity_dict[ftabgnps.loc[idx,'subclass_id']]
            elif ftabgnps.loc[idx, 'superclass_name']!='' and \
               ftabgnps.loc[idx, 'superclass_id']=='':
                ftabgnps.loc[idx, 'superclass_id'] = entity_dict[ftabgnps.loc[idx,'superclass_name']]
                if ftabgnps.loc[idx,'class_name']!='':
                    ftabgnps.loc[idx, 'class_id'] = entity_dict[ftabgnps.loc[idx,'class_name']]
                if ftabgnps.loc[idx,'subclass_name']!='':
                    ftabgnps.loc[idx, 'subclass_id'] = entity_dict[ftabgnps.loc[idx,'subclass_name']]
            else:
                pass
        ftabgnps.drop(['inchikey'], axis=1, inplace=True)
    else:
        ftabgnps.fillna('', inplace=True)
        ftabgnps['superclass_id'] = ''
        ftabgnps['class_id'] = ''
        ftabgnps['subclass_id'] = ''
        for idx in ftabgnps.index:
            if ftabgnps.loc[idx, 'superclass_name']!='':
                ftabgnps.loc[idx, 'superclass_id'] = entity_dict[ftabgnps.loc[idx,'superclass_name']]
                if ftabgnps.loc[idx,'class_name']!='':
                    ftabgnps.loc[idx, 'class_id'] = entity_dict[ftabgnps.loc[idx,'class_name']]
                if ftabgnps.loc[idx,'subclass_name']!='':
                    ftabgnps.loc[idx, 'subclass_id'] = entity_dict[ftabgnps.loc[idx,'subclass_name']]
            else:
                pass

    #inchikey = [Chem.InchiToInchiKey(x) for x in inchi]
    #class_anno = query_inchikey(inchikey)

    if outfile!='':
        ftabgnps.to_csv(os.path.join(doutput, outfile), index=None, sep='\t')
    return ftabgnps

def getLibraIDlink(i, tab):
    burl = 'https://gnps.ucsd.edu/ProteoSAFe/result.jsp?task='
    purl = '&view=view_all_annotations_DB#%7B%22main.Compound_Name_input%22%3A%22'
    eurl = '%22%7D'
    gid = tab.loc[0,'GNPSLinkout_Cluster'].split('&')[0].split('=')[1]
    if i!='':
        return '<a href=\"'+ burl+gid+purl+str(i)+eurl+ '\" target=\"_blank\">'+str(i)+'</a>'
    else:
        return ''

def getClusterSpectralink(i, tab):
    burl = 'https://gnps.ucsd.edu/ProteoSAFe/result.jsp?task='
    purl = '&view=view_all_clusters_withID&show=true#%7B%22main.cluster%20index_lowerinput%22%3A%22'
    eurl = '%22%2C%22main.cluster%20index_upperinput%22%3A%22'
    eeurl = '%22%7D'
    gid = tab.loc[0,'GNPSLinkout_Cluster'].split('&')[0].split('=')[1]
    if i!='':
        return '<a href=\"'+ burl+gid+purl+str(i)+eurl+ str(i)+eeurl+'\"target=\"_blank\">'+'N/A'+'</a>'
    else:
        return ''


def getLibraIDlink(i, tab):
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

