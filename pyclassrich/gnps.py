
import pandas as pd
import requests
import io

class Gnps:
    def __init__(self, taskid, workflow):
        """ Build GNPS node attributes and edge list
        Returns
        -------
        gnps  :  gnps task identification data
        """
        self.taskid = taskid
        self.workflow = workflow

    def description(self):
        return "{} taskid is a {} workflow".format(self.taskid, self.workflow)

    def getGnps(self): 
        """ Sends a request to ProteoSAFe.
        Parameters
        ----------
        taskid : str
           gnps task id
        workflow : str
           gnps workflow type
        Returns
        -------
        gnps : pandas.DataFrame
            node attributes table.
        net : pandas.DataFrame
            edge list
        """
        taskid = self.taskid
        taskid = taskid.split(',')
        workflow = self.workflow
        gdict = {}

        if workflow=='MZmine':
            url_to_attributes = "http://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task=%s&block=main&file=clusterinfo_summary/" % (taskid[0])
            url_to_edges = "http://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task=%s&block=main&file=networking_pairs_results_file_filtered/" % (taskid[0])
            url_to_features = "http://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task=%s&block=main&file=quantification_table/" % (taskid[0])
            gdict['gnps'] = pd.read_table(io.StringIO(requests.get(url_to_attributes).text))
            gdict['net'] = pd.read_table(io.StringIO(requests.get(url_to_edges).text))
            if len(taskid) > 1:
                url_to_attributes = "http://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task=%s&block=main&file=clusterinfo_summary/" % (taskid[1])
                url_to_edges = "http://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task=%s&block=main&file=networking_pairs_results_file_filtered/" % (taskid[1])
                gdict['gnps1'] = pd.read_table(io.StringIO(requests.get(url_to_attributes).text))
                gdict['net1'] = pd.read_table(io.StringIO(requests.get(url_to_edges).text))
            else:
                gdict['gnps1'] = None
                gdict['net1'] = None
        elif workflow=='V2':
            url_to_attributes = "http://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task=%s&block=main&file=clusterinfosummarygroup_attributes_withIDs_withcomponentID/" % (taskid[0])
            url_to_edges = "http://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task=%s&block=main&file=networkedges_selfloop/" % (taskid[0])
            gdict['gnps'] = pd.read_table(io.StringIO(requests.get(url_to_attributes).text))
            gdict['net'] = pd.read_table(io.StringIO(requests.get(url_to_edges).text))
            if len(taskid) > 1:
                url_to_attributes = "http://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task=%s&block=main&file=clusterinfosummarygroup_attributes_withIDs_withcomponentID/" % (taskid[1])
                url_to_edges = "http://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task=%s&block=main&file=networkedges_selfloop/" % (taskid[1])
                gdict['gnps1'] = pd.read_table(io.StringIO(requests.get(url_to_attributes).text))
                gdict['net1'] = pd.read_table(io.StringIO(requests.get(url_to_edges).text))
            else:
                gdict['gnps1'] = None
                gdict['net1'] = None

        return gdict
